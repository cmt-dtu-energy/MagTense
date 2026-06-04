from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from magtense.micromag import MicromagProblem

MU0 = 4 * np.pi * 1e-7
DEFAULT_MS = 1.61 / MU0  # [A/m], equivalent to mu0 Ms = 1.61 T.
DEFAULT_K0 = 4.3e6  # [J/m^3], uniaxial anisotropy constant.
DEFAULT_A0 = 7.7e-12  # [J/m], exchange stiffness used by the grain examples.
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
    t_end: float = 40e-9,
    nt: int = 2,
    fmm_cells_per_node: int = 660,
    fmm_eps: float = 1e-4,
    ifunif: int = 1,
    nlmin: int = 1,
    nlmax: int = 5,
    allow_fmm_short_circuit: int = 0,
    fmm_min_n: int = 20000,
    fmm_nterms: int = -1,
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
        solver="explicit",
        m0=np.tile(easy_axis, (ntot, 1)),
        A0=A0,
        Ms=Ms,
        K0=K0,
        alpha=4000.0,
        gamma=0.0,
        cuda=cuda,
        cvode=cvode,
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

    problem.log_dir = "timer_logs_single_grain"
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


def output_stem(L: float, n: int, use_fmm: bool) -> str:
    """Stable run identifier containing grain size, resolution and FMM state."""
    fmm_tag = "fmm_on" if use_fmm else "fmm_off"
    return f"single_grain_L{L:.6e}m_n{n}_{fmm_tag}"


def run_single_grain_coercivity(
    *,
    L: float,
    n: int,
    use_fmm: bool = False,
    cuda: bool = False,
    cvode: bool = False,
    tilt_degrees: float = DEFAULT_TILT_DEGREES,
    field_min_t: float = -7.0,
    field_max_t: float = 1.0,
    field_step_t: float = -0.1,
    output_dir: Path = Path("results"),
    plotting: bool = True,
    fmm_cells_per_node: int = 660,
    fmm_eps: float = 1e-4,
    ifunif: int = 1,
    nlmin: int = 1,
    nlmax: int = 5,
    allow_fmm_short_circuit: int = 0,
    fmm_min_n: int = 20000,
    fmm_nterms: int = -1,
) -> None:
    """Run one size/resolution hysteresis curve and save compatible arrays."""
    steps_t = np.arange(field_max_t, field_min_t + 0.5 * field_step_t, field_step_t)
    H_ext = build_hysteresis_field(steps_t)
    problem = create_single_grain_problem(
        L,
        n,
        use_fmm=use_fmm,
        cuda=cuda,
        cvode=cvode,
        tilt_degrees=tilt_degrees,
        fmm_cells_per_node=fmm_cells_per_node,
        fmm_eps=fmm_eps,
        ifunif=ifunif,
        nlmin=nlmin,
        nlmax=nlmax,
        allow_fmm_short_circuit=allow_fmm_short_circuit,
        fmm_min_n=fmm_min_n,
        fmm_nterms=fmm_nterms,
    )

    stem = output_stem(L, n, use_fmm)
    problem.timer_log_file = f"{stem}_timer.log"
    problem.trace_log_file = f"{stem}_trace.log"

    print("Single-grain coercivity experiment")
    print(f"  L = {L:.6e} m")
    print(f"  n = {n} ({n**3} cells)")
    print(f"  characteristic_length = {characteristic_length():.6e} m")
    print(f"  L / characteristic_length = {L / characteristic_length():.6g}")
    print(f"  easy_axis_tilt = {tilt_degrees:.3f} degrees")
    print(f"  easy_axis = {tilted_easy_axis(tilt_degrees)}")
    print(f"  use_fmm = {use_fmm}")
    print(f"  cuda = {cuda}")
    print(f"  cvode = {cvode}")

    start_time = time.time()
    res = problem.run_hysteresis(H_ext=H_ext)
    runtime = time.time() - start_time

    Mx, My, Mz = extract_mean_magnetisation(res, len(steps_t), DEFAULT_MS)
    H_A_per_m = H_ext[:, 3]
    H_T = MU0 * H_A_per_m
    Hc_A_per_m = interpolated_coercivity(H_A_per_m, Mz)
    Hc_T = MU0 * Hc_A_per_m if np.isfinite(Hc_A_per_m) else np.nan
    anisotropy_field_A_per_m = 2.0 * DEFAULT_K0 / DEFAULT_MS

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
        characteristic_length=characteristic_length(),
        size_factor=L / characteristic_length(),
        easy_axis=tilted_easy_axis(tilt_degrees),
        easy_axis_tilt_degrees=tilt_degrees,
        Ms=DEFAULT_MS,
        K0=DEFAULT_K0,
        A0=DEFAULT_A0,
    )
    print(f"Hysteresis simulation took {runtime:.3f} seconds")
    print(f"Hc = {Hc_A_per_m:.6e} A/m ({Hc_T:.6e} T)")
    print(f"Saved result object and scalars to: {res_path}")

    if plotting:
        plt = __import__("matplotlib.pyplot", fromlist=["pyplot"])

        fig_path = output_dir / f"{stem}.png"
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(H_A_per_m, Mz / DEFAULT_MS, ".-k", linewidth=1.5, markersize=4)
        ax.plot(Hc_A_per_m, 0.0, "r*", markersize=10, label="Hc")
        ax.set_xlabel("H_z [A/m]")
        ax.set_ylabel("M_z / Ms [-]")
        ax.set_title(f"Single grain hysteresis: L={L:.3e} m, n={n}, FMM={use_fmm}")
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.legend()
        fig.tight_layout()
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
    parser.add_argument("--use-fmm", action="store_true", help="Enable FMM demag.")
    parser.add_argument("--cuda", action="store_true", help="Enable CUDA acceleration.")
    parser.add_argument(
        "--cvode", action="store_true", help="Enable CVODE integration."
    )
    parser.add_argument(
        "--tilt-degrees",
        type=float,
        default=DEFAULT_TILT_DEGREES,
        help="Easy-axis tilt away from z in the x-z plane [degrees].",
    )
    parser.add_argument("--field-min-t", type=float, default=-7.0)
    parser.add_argument("--field-max-t", type=float, default=1.0)
    parser.add_argument("--field-step-t", type=float, default=-0.1)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable the default hysteresis plot written to the output directory.",
    )
    parser.add_argument("--fmm-cpn", type=int, default=660)
    parser.add_argument("--fmm-eps", type=float, default=1e-4)
    parser.add_argument("--ifunif", type=int, default=1)
    parser.add_argument("--nlmin", type=int, default=1)
    parser.add_argument("--nlmax", type=int, default=5)
    parser.add_argument("--allow-fmm-short-circuit", type=int, default=0)
    parser.add_argument("--fmm-min-n", type=int, default=20000)
    parser.add_argument("--fmm-nterms", type=int, default=-1)
    args = parser.parse_args()

    size_factor = 1.0 if args.size_factor is None else args.size_factor
    L = args.L if args.L is not None else size_factor * characteristic_length()
    run_single_grain_coercivity(
        L=L,
        n=args.n,
        use_fmm=args.use_fmm,
        cuda=args.cuda,
        cvode=args.cvode,
        tilt_degrees=args.tilt_degrees,
        field_min_t=args.field_min_t,
        field_max_t=args.field_max_t,
        field_step_t=args.field_step_t,
        output_dir=args.output_dir,
        plotting=not args.no_plot,
        fmm_cells_per_node=args.fmm_cpn,
        fmm_eps=args.fmm_eps,
        ifunif=args.ifunif,
        nlmin=args.nlmin,
        nlmax=args.nlmax,
        allow_fmm_short_circuit=args.allow_fmm_short_circuit,
        fmm_min_n=args.fmm_min_n,
        fmm_nterms=args.fmm_nterms,
    )


if __name__ == "__main__":
    main()
