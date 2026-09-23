"""Standard problem 6: domain wall pinning at a phase boundary.

Problem description
-------------------
A 1-D two-phase ferromagnet is modelled as a chain of ``x_steps`` cells along an 80 nm long
sample with a 1 nm x 1 nm cross section. The left half ("soft") and right half ("hard") differ
in exchange stiffness (A0), uniaxial anisotropy (K0) and / or saturation magnetisation (Ms) as
controlled by the ``settings`` string.

An antiparallel domain wall is initialised at the phase boundary. A field that increases
linearly with time is applied along the easy axis. At the depinning field the domain wall sweeps
through the hard phase and the full sample aligns with the field.

With ``use_minimizer`` the same field values are visited as a sequence of constant fields and the
equilibrium at each is found by the energy minimizer. That gives the static depinning field,
which is what the analytical values in the reference are, whereas the time ramp gives a
rate-dependent one.

The sample can be oriented along x, y or z with ``cart_dir``, which has no influence on the
result but checks the physics in every direction, and it can be meshed as a uniform grid or as a
chain of unstructured prisms with ``mesh_type``. ``two_d_sim`` extends the chain to a strip of
``two_d_size`` chains along y. This is the Python counterpart of
matlab/examples/Micromagnetism/mumag_micromag_Std_problem_6/Standard_problem_6.m and takes the
same options.

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

MESH_TYPES = ("uniform", "unstructuredPrisms")
SAMPLE_LENGTH = 80e-9  # m


def std_prob_6(
    settings: str = "akj",
    x_steps: int = 80,
    field_steps: int = 201,
    cart_dir: str = "x",
    mesh_type: str = "uniform",
    cuda: bool = True,
    cvode: bool = False,
    use_minimizer: bool = False,
    two_d_sim: bool = False,
    two_d_size: int = 5,
    plotting: bool = True,
    figpath: Path | None = None,
) -> float | None:
    """Simulate the muMag standard problem 6.

    Args:
        settings:    Which material parameters are reduced in the soft (left) half. Any subset
                     of ``'a'`` (exchange), ``'k'`` (anisotropy), ``'j'`` (Ms).
        x_steps:     Number of cells along the easy axis (default 80).
        field_steps: Number of time / field steps (default 201).
        cart_dir:    ``'x'``, ``'y'`` or ``'z'``, the axis the sample and the field lie along.
                     Has no influence on the result, but tests the physics in every direction.
        mesh_type:   ``'uniform'`` for a regular grid or ``'unstructuredPrisms'`` for the same
                     chain of cells given as unstructured prisms. The latter only works along x.
        cuda:        Enable GPU acceleration via CUDA.
        cvode:       Use CVODE time integrator instead of RK45.
        use_minimizer: Visit the same field values as constant fields and relax at each with the
                     energy minimizer instead of integrating the Landau-Lifshitz equation along
                     the time ramp. The number of effective-field evaluations spent is printed
                     either way.
        two_d_sim:   Run a 2D simulation: a strip of ``two_d_size`` chains along y, on the
                     uniform grid along x.
        two_d_size:  The number of chains of the 2D strip.
        plotting:    Show or save a plot of <m> vs applied field.
        figpath:     Directory for saving the figure. ``None`` shows it interactively.

    Returns:
        The depinning (switching) field in Tesla, or ``None`` if no switching occurred within the
        simulated field range.
    """
    if mesh_type not in MESH_TYPES:
        raise ValueError(f"mesh_type must be one of {MESH_TYPES}, got {mesh_type!r}")
    if cart_dir not in ("x", "y", "z"):
        raise ValueError(f"cart_dir must be 'x', 'y' or 'z', got {cart_dir!r}")
    if two_d_sim and (mesh_type != "uniform" or cart_dir != "x"):
        raise ValueError("the 2D simulation is only implemented on the uniform grid along x")
    if mesh_type == "unstructuredPrisms" and cart_dir != "x":
        raise ValueError("the unstructured mesh is only implemented along x")

    mu0 = 4 * np.pi * 1e-7
    dim = {"x": 0, "y": 1, "z": 2}[cart_dir]

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
    # The sample is 80 nm long along cart_dir, cut into x_steps cells, with a 1 nm cross section
    grid_L = np.array([1e-9, 1e-9, 1e-9])
    grid_L[dim] = SAMPLE_LENGTH
    n_half = x_steps // 2

    mesh_kwargs = {}
    if two_d_sim:
        res = (x_steps, two_d_size, 1)
        n_rows = two_d_size
    elif mesh_type == "uniform":
        res = [1, 1, 1]
        res[dim] = x_steps
        res = tuple(res)
        n_rows = 1
    else:
        # The same chain of cells, given as unstructured prisms
        res = (x_steps, 1, 1)
        n_rows = 1
        cell = SAMPLE_LENGTH / x_steps
        grid_pts = np.zeros((x_steps, 3))
        grid_pts[:, dim] = np.linspace(cell / 2, SAMPLE_LENGTH - cell / 2, x_steps)
        grid_abc = 1e-9 * np.ones((x_steps, 3))
        grid_abc[:, dim] = cell
        mesh_kwargs = {"grid_type": mesh_type, "grid_pts": grid_pts, "grid_abc": grid_abc}
    ntot = int(np.prod(res))

    # The cells of the soft (left) half. The tiles are ordered with x running fastest, so in the
    # 2D strip every chain is a block of x_steps consecutive tiles
    soft = np.concatenate([row * x_steps + np.arange(n_half) for row in range(n_rows)])

    # ── Per-cell material arrays: hard throughout, soft in left half ──────
    Ms_arr = Ms_hard * np.ones((ntot, 1))
    K0_arr = K0_hard * np.ones((ntot, 1))
    A0_arr = A0_hard * np.ones((ntot, 1))

    if "a" in settings:
        A0_arr[soft] = A0_soft
    if "k" in settings:
        K0_arr[soft] = K0_soft
    if "j" in settings:
        Ms_arr[soft] = Ms_soft

    # ── Problem setup ─────────────────────────────────────────────────────
    # The time ramp below is integrated by the dynamic solver. With the minimizer the same field
    # table is read as a sequence of constant fields and the equilibrium at each is found by the
    # energy minimizer, which gives the static depinning field - the quantity the analytical
    # values above are - rather than the rate-dependent one of the ramp.
    problem = MicromagProblem(
        res=res,
        grid_L=grid_L,
        solver="minimizer" if use_minimizer else "dynamic",
        alpha=alpha,
        gamma=gamma,
        Ms=Ms_arr,
        K0=K0_arr,
        A0=A0_arr,
        cuda=cuda,
        cvode=cvode,
        usedemag=False,   # demag negligible for this 1-D geometry
        usereturnhall=True,
        **mesh_kwargs,
    )
    problem.timer_log_file = "std_6_timer.log"
    problem.trace_log_file = "std_6_trace.log"

    # Easy axis along cart_dir for all cells
    easy_axis = np.zeros(3)
    easy_axis[dim] = 1.0
    problem.u_ea[:] = easy_axis

    # ── Initial state: antiparallel domain wall at the phase boundary ─────
    # Right half (hard): m along -cart_dir (anti-aligned with the field)
    # Left  half (soft): m along +cart_dir (aligned with the field)
    init_dir = {"x": [-1.0, 0.3, 0.0], "y": [0.0, -1.0, 0.3], "z": [0.0, 0.3, -1.0]}[cart_dir]
    init_dir = np.asarray(init_dir) / np.linalg.norm(init_dir)
    m0 = np.tile(init_dir, (ntot, 1))
    m0[soft, dim] = -m0[soft, dim]   # flip the soft region onto the field direction
    problem.m0 = m0

    # ── Applied field: linear ramp along +cart_dir ───────────────────────
    # μ₀·H(t) = 2×10⁷·(t + 20 ns) [T], i.e. 0.40 T at t = 0 and 2.40 T at t = 100 ns. The offset
    # is needed for the ill-conditioned 'k' and 'kj' cases and could be 0 in any other case.
    field_rate = 2e7 / mu0          # [A/m/s]
    H_offset = field_rate * 20e-9   # initial offset [A/m]

    def h_ext_fct(t: np.ndarray) -> np.ndarray:
        H = field_rate * np.asarray(t) + H_offset   # shape (nt_h_ext,)
        return np.outer(H, easy_axis)               # shape (nt_h_ext, 3)

    # ── Run ───────────────────────────────────────────────────────────────
    # All cells have the same volume on both meshes, so the mean over the cells is the mean
    # magnetization of the sample
    t_end = 100e-9
    t_fields = np.linspace(0, t_end, field_steps)
    if use_minimizer:
        # The field values the time ramp passes through at the output times, visited as a
        # sequence of constant fields. Column 0 is only a label; the solver reads the H vector.
        H_ext = np.zeros((field_steps, 4))
        H_ext[:, 0] = t_fields
        H_ext[:, 1:4] = h_ext_fct(t_fields)
        # The time window is only used if the minimizer has to fall back to the time integration
        problem.t = np.linspace(0, 10e-9, 2)
        problem.nt = 2
        M_out = problem.run_hysteresis(H_ext)[1]
        # The final state at each field, averaged over the cells
        M_avg = np.mean(M_out[-1, :, :, dim], axis=0)   # shape (field_steps,)
        mu0_H = mu0 * (field_rate * t_fields + H_offset)
        print(
            f'   settings="{settings}", minimizer: {int(problem.n_feval.sum())} field evaluations, '
            f"{int(problem.min_iter.sum())} iterations, "
            f"{int(np.sum(problem.min_status == 2))} fields not converged"
        )
    else:
        result = problem.run_simulation(
            t_end=t_end,
            nt=field_steps,
            fct_h_ext=h_ext_fct,
            nt_h_ext=field_steps,
        )
        t_out, M_out = result[0], result[1]

        # Average magnetisation along the easy axis at each output time
        M_sq = np.squeeze(M_out.copy(), axis=2)   # shape (nt, ntot, 3)
        M_avg = np.mean(M_sq[:, :, dim], axis=1)  # shape (nt,)

        # Applied field in Tesla at each output time
        mu0_H = mu0 * (field_rate * t_out + H_offset)
        print(f'   settings="{settings}", LL time ramp: {int(problem.n_feval.sum())} field evaluations')

    # Depinning field: the smallest field at which <m> > 1 − 10⁻³
    switched = M_avg > (1.0 - 1e-3)
    switching_field: float | None = float(mu0_H[switched].min()) if switched.any() else None

    if plotting:
        theory = THEORETICAL_PINNING_FIELDS.get(settings)
        _, ax = plt.subplots()
        ax.plot(mu0_H, M_avg, "b-", label=rf"$\langle m_{cart_dir} \rangle$")
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
        ax.set_ylabel(r"$\langle m \rangle$ [−]")
        method = "energy minimizer" if use_minimizer else "LL time ramp"
        ax.set_title(
            f'Standard problem 6, settings="{settings}", along {cart_dir}, mesh: {mesh_type}, {method}'
        )
        ax.legend()
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            suffix = "_minimizer" if use_minimizer else ""
            if cart_dir != "x":
                suffix += f"_{cart_dir}"
            if mesh_type != "uniform":
                suffix += f"_{mesh_type}"
            plt.savefig(figpath / f"6_settings_{settings}{suffix}.png")
        plt.close()

    return switching_field


if __name__ == "__main__":
    use_minimizer = False   # set True to relax with the energy minimizer at each field instead
    for _s in ("akj", "ak", "a", "k"):
        _sw = std_prob_6(
            settings=_s,
            x_steps=80,
            field_steps=201,
            cart_dir="x",
            mesh_type="uniform",   # 'uniform' or 'unstructuredPrisms'
            cuda=True,
            cvode=False,
            use_minimizer=use_minimizer,
            plotting=True,
            figpath=Path(__file__).resolve().parent,
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
