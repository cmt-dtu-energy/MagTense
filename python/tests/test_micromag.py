"""Tests for micromagnetic simulations using MagTense.

Standard problem 6
------------------
Domain wall pinning at a phase boundary.  A 1-D two-phase ferromagnet is
driven by a linearly increasing field; the test verifies that the depinning
(switching) field matches the analytical value to within 15 %.

Reference: Heistracher et al. (2022), https://arxiv.org/abs/2107.07855
"""

from __future__ import annotations

import numpy as np
import pytest

from magtense.micromag import MicromagProblem

# Analytical depinning fields [T] from Table I of the reference paper
_THEORETICAL_PINNING_FIELDS: dict[str, float] = {
    "akj": 1.568,
    "ak": 1.089,
    "aj": 1.206,
    "a": 0.838,
    "kj": 1.005,
    "k": 0.565,
    "j": 0.0,
    "": 0.0,
}


def _run_std_prob_6(
    settings: str = "a",
    x_steps: int = 80,
    field_steps: int = 201,
    cuda: bool = False,
    cvode: bool = False,
) -> float | None:
    """Run standard problem 6 and return the switching field in Tesla.

    Args:
        settings:    Subset of ``'a'`` (exchange), ``'k'`` (anisotropy),
                     ``'j'`` (Ms) indicating which parameters are reduced in
                     the soft (left) half of the sample.
        x_steps:     Number of cells along the easy axis.
        field_steps: Number of time / field steps.
        cuda:        Use GPU via CUDA.
        cvode:       Use CVODE time integrator.

    Returns:
        Depinning field in Tesla, or ``None`` if no switching occurred.
    """
    mu0 = 4 * np.pi * 1e-7

    # Hard-phase material parameters
    Ms_hard = 1.0 / mu0   # 1 T
    K0_hard = 1e6          # J/m³
    A0_hard = 1e-11        # J/m

    # Soft-phase material parameters
    Ms_soft = 0.25 / mu0
    K0_soft = 1e5
    A0_soft = 0.25e-11

    # Damping / gyromagnetic ratio (heavily damped: α₀ = 1)
    gamma0 = 2.2128e5
    alpha0 = 1.0
    gamma = gamma0 / (1 + alpha0 ** 2)
    alpha = alpha0 * gamma0 / (1 + alpha0 ** 2)

    n_half = x_steps // 2

    # Per-cell arrays: hard everywhere, then overwrite left half (soft region)
    Ms_arr = Ms_hard * np.ones((x_steps, 1))
    K0_arr = K0_hard * np.ones((x_steps, 1))
    A0_arr = A0_hard * np.ones((x_steps, 1))
    if "a" in settings:
        A0_arr[:n_half] = A0_soft
    if "k" in settings:
        K0_arr[:n_half] = K0_soft
    if "j" in settings:
        Ms_arr[:n_half] = Ms_soft

    problem = MicromagProblem(
        res=(x_steps, 1, 1),
        grid_L=[x_steps * 1e-9, 1e-9, 1e-9],
        alpha=alpha,
        gamma=gamma,
        Ms=Ms_arr,
        K0=K0_arr,
        A0=A0_arr,
        cuda=cuda,
        cvode=cvode,
        usedemag=False,   # demag negligible for this 1-D geometry
    )

    # Easy axis along x for all cells
    problem.u_ea[:, 0] = 1.0

    # Antiparallel domain wall at the phase boundary:
    #   right half (hard): mx ≈ −0.958  (anti-aligned with the +x field)
    #   left  half (soft): mx ≈ +0.958  (aligned with the +x field)
    init_dir = np.array([-1.0, 0.3, 0.0])
    init_dir /= np.linalg.norm(init_dir)
    m0 = np.tile(init_dir, (x_steps, 1))
    m0[:n_half, 0] = -m0[:n_half, 0]
    problem.m0 = m0

    # Linear field ramp along +x: μ₀·H(t) = 2×10⁷·(t + 20 ns) [T]
    field_rate = 2e7 / mu0           # [A/m/s]
    H_offset = field_rate * 20e-9    # initial offset [A/m]

    def h_ext_fct(t: np.ndarray) -> np.ndarray:
        H = field_rate * t + H_offset
        return np.outer(H, [1.0, 0.0, 0.0])

    t_end = 100e-9
    result = problem.run_simulation(
        t_end=t_end,
        nt=field_steps,
        fct_h_ext=h_ext_fct,
        nt_h_ext=field_steps,
    )
    t_out, M_out = result[0], result[1]

    # Average magnetisation along x at each output time
    M_sq = np.squeeze(M_out.copy(), axis=2)   # (nt, ntot, 3)
    M_avg = np.mean(M_sq[:, :, 0], axis=1)    # (nt,)

    # Applied field in Tesla at each output time
    mu0_H = 2e7 * (t_out + 20e-9)

    # Depinning field: first instant when ⟨m_x⟩ > 1 − 10⁻³
    switched = M_avg > (1.0 - 1e-3)
    return float(mu0_H[np.argmax(switched)]) if switched.any() else None


@pytest.mark.parametrize(
    ("settings", "theoretical"),
    [
        ("a", 0.838),    # exchange-only soft region; switches at ≈ 0.84 T
        ("akj", 1.568),  # fully soft region;         switches at ≈ 1.57 T
    ],
)
def test_std_prob_6(settings: str, theoretical: float) -> None:
    """Depinning field must be within 15 % of the analytical value."""
    switching_field = _run_std_prob_6(settings=settings, x_steps=80, field_steps=201)

    assert switching_field is not None, (
        f"No switching observed for settings='{settings}' "
        f"(expected ~{theoretical:.3f} T)"
    )

    rel_error = abs(switching_field - theoretical) / theoretical
    print(
        f"settings='{settings}': switching={switching_field:.4f} T, "
        f"theory={theoretical:.3f} T, rel_error={rel_error:.1%}"
    )
    assert rel_error < 0.15, (
        f"Depinning field {switching_field:.4f} T deviates {rel_error:.1%} "
        f"from theoretical {theoretical:.3f} T (tolerance 15 %)"
    )
