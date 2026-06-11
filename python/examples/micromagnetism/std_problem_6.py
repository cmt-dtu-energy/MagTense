"""Standard problem 6: domain wall pinning at a phase boundary.

Problem description
-------------------
A 1-D two-phase ferromagnet is modelled as a chain of ``x_steps`` cells,
each 1 nm × 1 nm × 1 nm.  The left half ("soft") and right half ("hard")
differ in exchange stiffness (A0), uniaxial anisotropy (K0) and / or
saturation magnetisation (Ms) as controlled by the ``settings`` string.

An antiparallel domain wall is initialised at the phase boundary.  A field
that increases linearly with time is applied along the easy axis.  At the
depinning field the domain wall sweeps through the hard phase and the full
sample aligns with the field.

Reference
---------
Heistracher et al., "Proposal for a micromagnetic standard problem:
domain wall pinning at phase boundaries" (2022).
https://arxiv.org/abs/2107.07855
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from magtense.micromag import MicromagProblem

# Analytical depinning fields [T] from the reference paper (Table I)
THEORETICAL_PINNING_FIELDS: dict[str, float] = {
    "akj": 1.568,
    "ak": 1.089,
    "aj": 1.206,
    "a": 0.838,
    "kj": 1.005,
    "k": 0.565,
    "j": 0.0,
    "": 0.0,
}


def std_prob_6(
    settings: str = "akj",
    x_steps: int = 80,
    field_steps: int = 201,
    cuda: bool = False,
    cvode: bool = False,
    plotting: bool = True,
    figpath: Path | None = None,
) -> float | None:
    """Simulate the muMag standard problem 6.

    Args:
        settings:    Which material parameters are reduced in the soft (left)
                     half.  Any subset of ``'a'`` (exchange), ``'k'``
                     (anisotropy), ``'j'`` (Ms).
        x_steps:     Number of cells along the easy axis (default 80).
        field_steps: Number of time / field steps (default 201).
        cuda:        Enable GPU acceleration via CUDA.
        cvode:       Use CVODE time integrator instead of RK45.
        plotting:    Show or save a plot of ⟨m_x⟩ vs applied field.
        figpath:     Directory for saving the figure.  ``None`` → show
                     interactively.

    Returns:
        The depinning (switching) field in Tesla, or ``None`` if no
        switching occurred within the simulated field range.
    """
    mu0 = 4 * np.pi * 1e-7

    # ── Material parameters ───────────────────────────────────────────────
    Ms_hard = 1.0 / mu0      # ≈ 795 775 A/m  (1 T)
    K0_hard = 1e6             # J/m³
    A0_hard = 1e-11           # J/m

    Ms_soft = 0.25 / mu0     # 0.25 T
    K0_soft = 1e5             # J/m³
    A0_soft = 0.25e-11        # J/m

    # Damping / gyromagnetic ratio (α₀ = 1, heavily damped)
    gamma0 = 2.2128e5
    alpha0 = 1.0
    gamma = gamma0 / (1 + alpha0 ** 2)
    alpha = alpha0 * gamma0 / (1 + alpha0 ** 2)

    # ── Geometry ──────────────────────────────────────────────────────────
    res = (x_steps, 1, 1)
    grid_L = [x_steps * 1e-9, 1e-9, 1e-9]
    n_half = x_steps // 2

    # ── Per-cell material arrays: hard throughout, soft in left half ──────
    Ms_arr = Ms_hard * np.ones((x_steps, 1))
    K0_arr = K0_hard * np.ones((x_steps, 1))
    A0_arr = A0_hard * np.ones((x_steps, 1))

    if "a" in settings:
        A0_arr[:n_half] = A0_soft
    if "k" in settings:
        K0_arr[:n_half] = K0_soft
    if "j" in settings:
        Ms_arr[:n_half] = Ms_soft

    # ── Problem setup ─────────────────────────────────────────────────────
    problem = MicromagProblem(
        res=res,
        grid_L=grid_L,
        alpha=alpha,
        gamma=gamma,
        Ms=Ms_arr,
        K0=K0_arr,
        A0=A0_arr,
        cuda=cuda,
        cvode=cvode,
        usedemag=False,   # demag negligible for this 1-D geometry
    )




    #--------- disable fmm -----
    problem.use_fmm = 0
    #----------------------------
    
    #--------------- set trace/timing options -------------
    problem.use_fmm = 0
    problem.window_enabled = 0
    problem.window_interval = 30.0
    problem.trace_enabled = 0
    problem.flush_each = 1
    problem.trace_verbose = 2
    problem.timer_log_file =  "std_6_timer.log"
    problem.trace_log_file =  "std_6_trace.log"
    #------------------------------------------------------



    # Easy axis along x for all cells
    problem.u_ea[:, 0] = 1.0

    # ── Initial state: antiparallel domain wall at the phase boundary ─────
    # Right half (hard): mx ≈ −0.958  (anti-aligned with the +x field)
    # Left  half (soft): mx ≈ +0.958  (aligned with the +x field)
    init_dir = np.array([-1.0, 0.3, 0.0])
    init_dir /= np.linalg.norm(init_dir)
    m0 = np.tile(init_dir, (x_steps, 1))
    m0[:n_half, 0] = -m0[:n_half, 0]   # flip soft-region x-component to +x
    problem.m0 = m0

    # ── Applied field: linear ramp along +x ──────────────────────────────
    # μ₀·H(t) = 2×10⁷·(t + 20 ns) [T]
    # At t = 0:     μ₀H = 0.40 T
    # At t = 100 ns: μ₀H = 2.40 T
    field_rate = 2e7 / mu0          # [A/m/s]
    H_offset = field_rate * 20e-9   # initial offset [A/m]

    def h_ext_fct(t: np.ndarray) -> np.ndarray:
        #H = field_rate * t + H_offset      # shape (nt_h_ext,)
        H = 2e7/mu0 * t + 2e7/mu0 * 20e-9  # shape (nt_h_ext,)
        return np.outer(H, [1.0, 0.0, 0.0])  # shape (nt_h_ext, 3)

    # ── Run ───────────────────────────────────────────────────────────────
    t_end = 100e-9
    result = problem.run_simulation(
        t_end=t_end,
        nt=field_steps,
        fct_h_ext=h_ext_fct,
        nt_h_ext=2001,
    )
    t_out, M_out = result[0], result[1]

    # ── Post-processing ───────────────────────────────────────────────────
    # Average magnetisation along x (easy axis) at each output time
    M_sq = np.squeeze(M_out.copy(), axis=2)   # shape (nt, ntot, 3)
    M_avg = np.mean(M_sq[:, :, 0], axis=1)    # shape (nt,)

    # Applied field in Tesla at each output time
    mu0_H = 4*np.pi*1e-7*(field_rate * t_out + H_offset)

    # Depinning field: first instant when ⟨m_x⟩ > 1 − 10⁻³
    switched = M_avg > (1.0 - 1e-3)
    switching_field: float | None = (
        float(mu0_H[np.argmax(switched)]) if switched.any() else None
    )

    if plotting:
        theory = THEORETICAL_PINNING_FIELDS.get(settings)
        _, ax = plt.subplots()
        ax.plot(mu0_H, M_avg, "b-", label=r"$\langle m_x \rangle$")
        if switching_field is not None:
            ax.axvline(
                switching_field,
                color="r",
                linestyle="--",
                label=f"Switching field = {switching_field:.3f} T",
            )
        if theory is not None:
            ax.axvline(
                theory,
                color="g",
                linestyle=":",
                label=f"Theoretical = {theory:.3f} T",
            )
        ax.set_xlabel(r"$\mu_0 H_{\rm app}$ [T]")
        ax.set_ylabel(r"$\langle m_x \rangle$ [−]")
        ax.set_title(f'Standard problem 6, settings="{settings}"')
        ax.legend()
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            plt.savefig(figpath / f"6_settings_{settings}.png")
        plt.close()

    return switching_field


if __name__ == "__main__":
    for _s in ("akj", "ak", "a", "k"):
    # iterate over the single settings string "akj" (not its characters)
    #for _s in ("akj",):
        _sw = std_prob_6(
            settings=_s,
            x_steps=80,
            field_steps=201,
            cuda=False,
            cvode=False,
            plotting=True,
            #figpath=None,
            figpath= Path.cwd() / "results"
        )
        _theory = THEORETICAL_PINNING_FIELDS[_s]
        if _sw is not None:
            print(
                f'settings="{_s}": switching field = {_sw:.4f} T  '
                f"(theory {_theory:.3f} T, "
                f"error {abs(_sw - _theory) / _theory * 100:.1f} %)"
            )
        else:
            print(f'settings="{_s}": no switching observed (theory {_theory:.3f} T)')
