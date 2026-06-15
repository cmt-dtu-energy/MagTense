from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from magtense.micromag import MicromagProblem

MU0 = 4 * np.pi * 1e-7
# DEFAULT_MS = 1.61 / MU0  # [A/m], equivalent to mu0 Ms = 1.61 T.
# DEFAULT_K0 = 4.3e6  # [J/m^3], uniaxial anisotropy constant.
# DEFAULT_A0 = 7.7e-12  # [J/m], exchange stiffness used by the grain examples.

#DEFAULT_MS = 2.3e6      # A/m
#DEFAULT_A0 = 7.0e-12    # J/m
#DEFAULT_K0 = 1.8e6      # J/m^3

#DEFAULT_MS = 1.25e6     # A/m


DEFAULT_MU0_MS_T = 2.4  # [T]
DEFAULT_MS = DEFAULT_MU0_MS_T / MU0  # [A/m]
DEFAULT_K0 = 1.0e6  # [J/m^3]
DEFAULT_A0 = 7.0e-12  # [J/m]

DEFAULT_TILT_DEGREES = 3.0


def characteristic_length(A0: float = DEFAULT_A0, Ms: float = DEFAULT_MS) -> float:
    """Return sqrt(A0 / (0.5 * mu0 * Ms**2)) in SI units [m]."""
    return float(np.sqrt(A0 / (0.5 * MU0 * Ms**2)))


def tilted_easy_axis(tilt_degrees: float = DEFAULT_TILT_DEGREES) -> np.ndarray:
    """
    Construct a unit easy-axis vector tilted in the x-z plane.

    The external hysteresis field is applied along +z.  A small default tilt of
    three degrees avoids the perfectly aligned switching geometry while keeping
    this a deliberately simple single-grain experiment.
    """
    theta = np.deg2rad(tilt_degrees)
    return np.array([np.sin(theta), 0.0, np.cos(theta)], dtype=float)


def build_hysteresis_field(steps_t: np.ndarray) -> np.ndarray:
    """
    Build MagTense hysteresis input with the field along z.

    The first column stores the signed scalar sweep value.  The vector columns
    are H in A/m, so the Tesla-equivalent sweep values are divided by mu0, which
    is the convention used in the existing grain hysteresis scripts.
    """
    h_ext = np.zeros((len(steps_t), 4), dtype=float)
    h_ext[:, 0] = steps_t
    h_ext[:, 3] = steps_t / MU0
    return h_ext


def create_single_grain_problem(
    L: float,
    n: int,
    *,
    use_fmm: bool = False,
    cuda: bool = False,
    cvode: bool = False,
    tilt_degrees: float = DEFAULT_TILT_DEGREES,
    Ms: float = DEFAULT_MS,
    K0: float = DEFAULT_K0,
    A0: float = DEFAULT_A0,
    t_end: float = 1e-9,
    nt: int = 2,
    fmm_cells_per_node: int = 660,
    fmm_eps: float = 1e-4,
    ifunif: int = 1,
    nlmin: int = 1,
    nlmax: int = 5,
    allow_fmm_short_circuit: int = 0,
    fmm_min_n: int = 20000,
    fmm_nterms: int = -1,
    timer_log_dir: Path = Path("timer_logs_single_grain"),
    hysteresis_solver: str = "static",
) -> MicromagProblem:
    """
    Create the single-grain cubic coercivity problem.

    Physical setup
    --------------
    * Uniform cubic domain with Lx = Ly = Lz = L [m].
    * Cubic numerical grid with nx = ny = nz = n.
    * All cells are one material/grain: identical A0, K0, Ms and easy axis.
    * External field is applied along z.
    * Easy axis is tilted slightly away from z in the x-z plane.

    The characteristic exchange/demagnetisation scale used by the run script is
    sqrt(A0 / (0.5 * mu0 * Ms**2)); default size factors are centred around this
    value and are passed as L = factor * characteristic_length.
    """
    ntot = n**3
    easy_axis = tilted_easy_axis(tilt_degrees)

    problem = MicromagProblem(
        res=[n, n, n],
        grid_L=[L, L, L],
        grid_type="uniform",
        grid_pts = None,
        grid_abc = None,
        solver="explicit",
        hysteresis_solver=hysteresis_solver,
        m0=np.tile([0.0, 0.0, 1.0], (ntot, 1)),
        #m0=np.tile(easy_axis, (ntot, 1)),
        A0=A0,
        Ms=Ms,
        K0=K0,
        alpha=4000.0,
        gamma=0.0,
        cuda=cuda,
        cvode=cvode,
        useavgn=1
    )

    problem.u_ea[:, :] = easy_axis[np.newaxis, :]
    problem.t = np.linspace(0.0, t_end, nt)
    problem.nt = len(problem.t)
    problem.t_conv = np.linspace(0.0, t_end, nt)
    problem.nt_conv = len(problem.t_conv)

    # A uniform cubic cell has up to 26 neighbours in the exchange stencil.
    # Keep this explicit so the experiment remains easy to extend later.
    problem.exch_presize = problem.ntot * (27 + 1)

    # FMM is intentionally disabled by default for this experiment.  The CLI and
    # run script can enable it later without changing the physical setup.
    problem.use_fmm = int(use_fmm)
    problem.fmm_cells_per_node = fmm_cells_per_node if use_fmm else 0
    problem.fmm_eps = fmm_eps
    problem.ifunif = ifunif
    problem.nlmin = nlmin
    problem.nlmax = nlmax
    problem.allow_fmm_short_circuit = allow_fmm_short_circuit
    problem.fmm_min_n = fmm_min_n
    problem.fmm_nterms = fmm_nterms

    problem.log_dir = str(timer_log_dir)
    problem.window_enabled = 1
    problem.window_interval = 30.0
    problem.trace_enabled = 0
    problem.flush_each = 1
    problem.trace_verbose = 1

    return problem


def extract_mean_magnetisation(
    res: list[np.ndarray | int],
    n_steps: int,
    Ms: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return volume-averaged magnetisation components in A/m."""
    M_out = res[1][1, :, :, :]  # Existing grain scripts use this hysteresis slice.
    Mx = Ms * np.mean(M_out[:, :n_steps, 0], axis=0)
    My = Ms * np.mean(M_out[:, :n_steps, 1], axis=0)
    Mz = Ms * np.mean(M_out[:, :n_steps, 2], axis=0)
    return Mx, My, Mz


def interpolated_coercivity(H: np.ndarray, M: np.ndarray) -> float:
    """Interpolate H(M=0) for the descending hysteresis branch."""
    H = np.asarray(H, dtype=float).ravel()
    M = np.asarray(M, dtype=float).ravel()
    crossings = np.flatnonzero(M[:-1] * M[1:] <= 0.0)
    if crossings.size == 0:
        return float("nan")

    i = int(crossings[0])
    m0, m1 = M[i], M[i + 1]
    h0, h1 = H[i], H[i + 1]
    if m1 == m0:
        return float(0.5 * (h0 + h1))
    return float(h0 + (0.0 - m0) * (h1 - h0) / (m1 - m0))


def output_stem(
    size_factor: float,
    n: int,
    use_fmm: bool,
    fmm_nterms: int,
    nlmax: int,
) -> str:
    """Stable run identifier containing size factor, resolution and FMM state."""
    stem = f"single_grain_sf{size_factor:.2g}_n{n}"
    if use_fmm:
        stem += f"_fmm_on_N{fmm_nterms}_L{nlmax}"
    return stem


def run_single_grain_coercivity(
    *,
    L: float,
    size_factor: float,
    n: int,
    use_fmm: bool = False,
    periodic: bool = False,
    cuda: bool = False,
    cvode: bool = False,
    tilt_degrees: float = DEFAULT_TILT_DEGREES,
    field_min_t: float = -7.0,
    field_max_t: float = 1.0,
    field_step_t: float = -0.1,
    output_dir: Path = Path("results"),
    timer_log_dir: Path = Path("timer_logs_single_grain"),
    plotting: bool = True,
    Ms: float = DEFAULT_MS,
    K0: float = DEFAULT_K0,
    A0: float = DEFAULT_A0,
    fmm_cells_per_node: int = 660,
    fmm_eps: float = 1e-4,
    ifunif: int = 1,
    nlmin: int = 1,
    nlmax: int = 5,
    allow_fmm_short_circuit: int = 0,
    fmm_min_n: int = 20000,
    fmm_nterms: int = -1,
    adaptive: bool = False,
    adaptive_max_steps: int | None = None,
    adaptive_dh_initial_t: float | None = None,
    adaptive_dh_min_t: float | None = None,
    adaptive_dh_max_t: float | None = None,
) -> None:
    """Run one size/resolution hysteresis curve and save compatible arrays."""
    timer_log_dir.mkdir(parents=True, exist_ok=True)
    steps_t = np.arange(field_max_t, field_min_t + 0.5 * field_step_t, field_step_t)
    H_ext = build_hysteresis_field(steps_t)
    problem = create_single_grain_problem(
        L,
        n,
        use_fmm=use_fmm,
        cuda=cuda,
        cvode=cvode,
        tilt_degrees=tilt_degrees,
        Ms=Ms,
        K0=K0,
        A0=A0,
        fmm_cells_per_node=fmm_cells_per_node,
        fmm_eps=fmm_eps,
        ifunif=ifunif,
        nlmin=nlmin,
        nlmax=nlmax,
        allow_fmm_short_circuit=allow_fmm_short_circuit,
        fmm_min_n=fmm_min_n,
        fmm_nterms=fmm_nterms,
        timer_log_dir=timer_log_dir,
        hysteresis_solver="adaptive" if adaptive else "static",
    )

    # Configure periodic exchange boundary conditions at the script level
    # (do not modify core MicromagProblem implementation here).
    if periodic:
        try:
            problem.exchPBC = [1, 1, 1]
        except Exception:
            # Best-effort: if attribute not present, set attribute anyway
            setattr(problem, "exchPBC", [1, 1, 1])
    else:
        try:
            problem.exchPBC = [0, 0, 0]
        except Exception:
            setattr(problem, "exchPBC", [0, 0, 0])

    # Also set periodic repetition for the demag field (n_macro):
    # use a single periodic repetition in each dimension when periodic=True.
    if periodic:
        try:
            problem.n_macro = np.array([1, 1, 1], dtype=int)
            problem.shiftVec = np.array([L, L, L])
        except Exception:
            setattr(problem, "n_macro", np.array([1, 1, 1], dtype=int))
            setattr(problem, "shiftVec", np.array([L, L, L]))
    else:
        try:
            problem.n_macro = np.zeros(3, dtype=int)
            problem.shiftVec = np.zeros(3)
        except Exception:
            setattr(problem, "n_macro", np.zeros(3, dtype=int))
            setattr(problem, "shiftVec", np.zeros(3))

    stem = output_stem(size_factor, n, use_fmm, fmm_nterms, nlmax)
    # Append _P to the stem when periodic boundary conditions are enabled.
    if periodic:
        n_token = f"_n{n}"
        if n_token in stem:
            stem = stem.replace(n_token, n_token + "_P", 1)
        else:
            stem = stem + "_P"
    if adaptive:
        dH_min_t = (
            abs(adaptive_dh_min_t)
            if adaptive_dh_min_t is not None
            else abs(field_step_t / 10.0)
        )
        stem += f"_A_FS{dH_min_t:.1e}"
    problem.timer_log_file = f"{stem}_timer.log"
    problem.trace_log_file = f"{stem}_trace.log"

    print("Single-grain coercivity experiment")
    print(f"  L = {L:.6e} m")
    print(f"  size_factor = {size_factor:.2g}")
    print(f"  n = {n} ({n**3} cells)")
    material_length = characteristic_length(A0, Ms)
    print(f"  mu0 Ms = {MU0 * Ms:.6e} T")
    print(f"  Ms = {Ms:.6e} A/m")
    print(f"  A0 = {A0:.6e} J/m")
    print(f"  K0 = {K0:.6e} J/m^3")
    print(f"  characteristic_length = {material_length:.6e} m")
    print(f"  L / characteristic_length = {L / material_length:.6g}")
    print(f"  easy_axis_tilt = {tilt_degrees:.3f} degrees")
    print(f"  easy_axis = {tilted_easy_axis(tilt_degrees)}")
    print(f"  use_fmm = {use_fmm}")
    print(f"  cuda = {cuda}")
    print(f"  cvode = {cvode}")
    print(f"  adaptive = {adaptive}")
    start_time = time.time()
    if adaptive:
        max_steps = (
            adaptive_max_steps if adaptive_max_steps is not None else len(steps_t)
        )
        dH_initial_t = (
            adaptive_dh_initial_t
            if adaptive_dh_initial_t is not None
            else field_step_t
        )
        dH_min_t = (
            adaptive_dh_min_t
            if adaptive_dh_min_t is not None
            else field_step_t / 10.0
        )
        dH_max_t = (
            adaptive_dh_max_t
            if adaptive_dh_max_t is not None
            else field_step_t * 2.0
        )
        dH_initial = abs(dH_initial_t) / MU0
        dH_min = abs(dH_min_t) / MU0
        dH_max = abs(dH_max_t) / MU0
        res = problem.run_hysteresis_adaptive(
            H_start=H_ext[0, 1:4],
            H_end=H_ext[-1, 1:4],
            dH_initial=dH_initial,
            dH_min=dH_min,
            dH_max=dH_max,
            max_steps=max_steps,
            switch_refine_dH=dH_min,
        )
        n_fields = res[-1]
        H_A_per_m = res[4][0, 0, :n_fields, 2]
    else:
        res = problem.run_hysteresis(H_ext=H_ext)
        n_fields = len(steps_t)
        H_A_per_m = H_ext[:, 3]
    runtime = time.time() - start_time

    Mx, My, Mz = extract_mean_magnetisation(res, n_fields, Ms)
    H_T = MU0 * H_A_per_m
    Hc_A_per_m = interpolated_coercivity(H_A_per_m, Mz)
    Hc_T = MU0 * Hc_A_per_m if np.isfinite(Hc_A_per_m) else np.nan
    anisotropy_field_A_per_m = 2.0 * K0 / Ms

    output_dir.mkdir(parents=True, exist_ok=True)
    res_path = output_dir / f"{stem}.npz"
    np.savez(
        res_path,
        res=np.array(res, dtype=object),
        H_array=H_T,
        H_array_A_per_m=H_A_per_m,
        M_array=Mz,
        Mx_array=Mx,
        My_array=My,
        Mz_array=Mz,
        Hc=Hc_A_per_m,
        Hc_A_per_m=Hc_A_per_m,
        Hc_T=Hc_T,
        H_N=anisotropy_field_A_per_m,
        runtime=runtime,
        L=L,
        n=n,
        ntot=n**3,
        use_fmm=use_fmm,
        cuda=cuda,
        cvode=cvode,
        fmm_cells_per_node=fmm_cells_per_node,
        fmm_eps=fmm_eps,
        ifunif=ifunif,
        nlmin=nlmin,
        nlmax=nlmax,
        allow_fmm_short_circuit=allow_fmm_short_circuit,
        fmm_min_n=fmm_min_n,
        fmm_nterms=fmm_nterms,
        characteristic_length=material_length,
        size_factor=L / material_length,
        easy_axis=tilted_easy_axis(tilt_degrees),
        easy_axis_tilt_degrees=tilt_degrees,
        mu0_Ms_T=MU0 * Ms,
        Ms=Ms,
        K0=K0,
        A0=A0,
    )
    print(f"Hysteresis simulation took {runtime:.3f} seconds")
    print(f"Hc = {Hc_A_per_m:.6e} A/m ({Hc_T:.6e} T)")
    print(f"Saved result object and scalars to: {res_path}")

    if plotting:
        plt = __import__("matplotlib.pyplot", fromlist=["pyplot"])

        fig_path = output_dir / f"{stem}.png"
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(MU0 * H_A_per_m, MU0 * Mz, ".-k", linewidth=1.5, markersize=4)
        ax.plot(
            MU0 * Hc_A_per_m,
            0.0,
            "r*",
            markersize=10,
            label=f"Hc = {Hc_T:.3f} T",
        )
        ax.set_xlabel(r"Applied Field, $\mu_0 H$ [ T ]")
        ax.set_ylabel(r"Magnetization, $\mu_0 M$ [ T ]")
        ax.set_title(
            f"Single grain hysteresis: size_factor={size_factor:.2g}, n={n}, FMM={use_fmm}"
        )
        ax.grid(True, which="both", linestyle=":", linewidth=0.8, alpha=0.6)
        ax.legend(loc="lower right", frameon=True, framealpha=0.9, edgecolor="inherit")
        fig.tight_layout(pad=0.2)
        fig.savefig(fig_path, dpi=300)
        plt.close(fig)
        print(f"Saved figure to: {fig_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run a simple single-grain grain-size/resolution coercivity test. "
            "With no arguments, this runs size_factor=1, n=10, FMM/CUDA/CVODE off, "
            "and writes the result plus plot to ./results."
        )
    )
    size_group = parser.add_mutually_exclusive_group()
    size_group.add_argument("--L", type=float, help="Cubic domain side length [m].")
    size_group.add_argument(
        "--size-factor",
        type=float,
        help=(
            "Set L = size_factor * sqrt(A0 / (0.5 * mu0 * Ms**2)); "
            "defaults to 1.0 when neither --L nor --size-factor is provided."
        ),
    )
    parser.add_argument(
        "--n", type=int, default=10, help="Cubic grid resolution (default: 10)."
    )
    parser.add_argument(
        "--tilt-degrees",
        type=float,
        default=DEFAULT_TILT_DEGREES,
        help="Easy-axis tilt away from z in the x-z plane [degrees].",
    )
    parser.add_argument("--field-min-t", type=float, default=-3.0)
    parser.add_argument("--field-max-t", type=float, default=1.0)
    parser.add_argument("--field-step-t", type=float, default=-0.1)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--timer-log-dir",
        type=Path,
        default=Path("timer_logs_single_grain"),
        help="Directory for backend timer and trace logs.",
    )
    parser.add_argument(
        "--mu0-ms-t",
        type=float,
        default=DEFAULT_MU0_MS_T,
        help="Saturation magnetization expressed as mu0*Ms [T].",
    )
    parser.add_argument(
        "--a0", type=float, default=DEFAULT_A0, help="Exchange stiffness [J/m]."
    )
    parser.add_argument(
        "--k0",
        type=float,
        default=DEFAULT_K0,
        help="Uniaxial anisotropy constant [J/m^3].",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable the default hysteresis plot written to the output directory.",
    )
    parser.add_argument("--use-fmm", action="store_true", help="Enable FMM demag.")
    parser.add_argument("--fmm-cpn", type=int, default=660)
    parser.add_argument("--fmm-eps", type=float, default=1e-4)
    parser.add_argument("--ifunif", type=int, default=1)
    parser.add_argument("--nlmin", type=int, default=1)
    parser.add_argument("--nlmax", type=int, default=5)
    parser.add_argument("--allow-fmm-short-circuit", type=int, default=0)
    parser.add_argument("--fmm-min-n", type=int, default=20000)
    parser.add_argument("--fmm-nterms", type=int, default=-1)
    parser.add_argument(
        "--adaptive",
        action="store_true",
        help="Use adaptive backend hysteresis stepping.",
    )
    parser.add_argument("--adaptive-max-steps", type=int, default=500)
    parser.add_argument("--adaptive-dh-initial-t", type=float, default=0.5)
    parser.add_argument("--adaptive-dh-min-t", type=float, default=0.01)
    parser.add_argument("--adaptive-dh-max-t", type=float, default=0.5)
    parser.add_argument(
        "--periodic",
        action="store_true",
        help="Enable periodic exchange boundary conditions (exchPBC=1) and append _P to output names.",
    )
    args = parser.parse_args()

    if not np.isfinite(args.mu0_ms_t) or args.mu0_ms_t <= 0.0:
        parser.error("--mu0-ms-t must be positive")
    if not np.isfinite(args.a0) or args.a0 <= 0.0:
        parser.error("--a0 must be positive")
    if not np.isfinite(args.k0) or args.k0 < 0.0:
        parser.error("--k0 must be non-negative")
    if args.n <= 0:
        parser.error("--n must be positive")
    if args.size_factor is not None and args.size_factor <= 0.0:
        parser.error("--size-factor must be positive")
    if args.L is not None and args.L <= 0.0:
        parser.error("--L must be positive")

    Ms = args.mu0_ms_t / MU0
    material_length = characteristic_length(args.a0, Ms)

    if args.size_factor is None:
        if args.L is None:
            size_factor = 1.0
        else:
            size_factor = args.L / material_length
    else:
        size_factor = args.size_factor
    L = args.L if args.L is not None else size_factor * material_length
    run_single_grain_coercivity(
        L=L,
        size_factor=size_factor,
        n=args.n,
        use_fmm=args.use_fmm,
        cuda=True,
        cvode=False,
        tilt_degrees=args.tilt_degrees,
        field_min_t=args.field_min_t,
        field_max_t=args.field_max_t,
        field_step_t=args.field_step_t,
        output_dir=args.output_dir,
        timer_log_dir=args.timer_log_dir,
        plotting=not args.no_plot,
        Ms=Ms,
        K0=args.k0,
        A0=args.a0,
        periodic=args.periodic,
        fmm_cells_per_node=args.fmm_cpn,
        fmm_eps=args.fmm_eps,
        ifunif=args.ifunif,
        nlmin=args.nlmin,
        nlmax=args.nlmax,
        allow_fmm_short_circuit=args.allow_fmm_short_circuit,
        fmm_min_n=args.fmm_min_n,
        fmm_nterms=args.fmm_nterms,
        adaptive=args.adaptive,
        adaptive_max_steps=args.adaptive_max_steps,
        adaptive_dh_initial_t=args.adaptive_dh_initial_t,
        adaptive_dh_min_t=args.adaptive_dh_min_t,
        adaptive_dh_max_t=args.adaptive_dh_max_t,
    )


if __name__ == "__main__":
    main()
