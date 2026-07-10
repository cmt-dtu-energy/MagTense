"""Construction and low-level processing helpers for single-grain runs."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from magtense.micromag import MicromagProblem

from .metrics import MU0, interpolated_coercivity

DEFAULT_MU0_MS_T = 2.4
DEFAULT_MS = DEFAULT_MU0_MS_T / MU0
DEFAULT_K0 = 1.0e6
DEFAULT_A0 = 7.0e-12
DEFAULT_TILT_DEGREES = 3.0


def characteristic_length(A0: float = DEFAULT_A0, Ms: float = DEFAULT_MS) -> float:
    """Return sqrt(A0 / (0.5 * mu0 * Ms**2)) in metres."""
    return float(np.sqrt(A0 / (0.5 * MU0 * Ms**2)))


def tilted_easy_axis(tilt_degrees: float = DEFAULT_TILT_DEGREES) -> np.ndarray:
    """Construct a unit easy-axis vector tilted in the x-z plane."""
    theta = np.deg2rad(tilt_degrees)
    return np.array([np.sin(theta), 0.0, np.cos(theta)], dtype=float)


def z_easy_axis() -> np.ndarray:
    """Return the fixed uniaxial easy axis used by the single-grain runs."""
    return np.array([0.0, 0.0, 1.0], dtype=float)


def tilted_field_direction(tilt_degrees: float = DEFAULT_TILT_DEGREES) -> np.ndarray:
    """Construct a unit external-field direction tilted in the x-z plane."""
    return tilted_easy_axis(tilt_degrees)


def build_hysteresis_field(
    steps_t: np.ndarray,
    direction: np.ndarray | None = None,
) -> np.ndarray:
    """Build the MagTense four-column external-field input."""
    direction = z_easy_axis() if direction is None else np.asarray(direction, dtype=float)
    norm = np.linalg.norm(direction)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("field direction must be a finite non-zero vector")
    direction = direction / norm
    h_ext = np.zeros((len(steps_t), 4), dtype=float)
    h_ext[:, 0] = steps_t
    h_ext[:, 1:4] = (steps_t / MU0)[:, np.newaxis] * direction[np.newaxis, :]
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
    easy_axis: np.ndarray | None = None,
    m0_direction: np.ndarray | None = None,
) -> MicromagProblem:
    """Create the uniform cubic micromagnetic problem used by this experiment."""
    ntot = n**3
    easy_axis = z_easy_axis() if easy_axis is None else np.asarray(easy_axis, dtype=float)
    m0_direction = (
        tilted_field_direction(tilt_degrees)
        if m0_direction is None
        else np.asarray(m0_direction, dtype=float)
    )
    problem = MicromagProblem(
        res=[n, n, n],
        grid_L=[L, L, L],
        grid_type="uniform",
        grid_pts=None,
        grid_abc=None,
        solver="explicit",
        hysteresis_solver=hysteresis_solver,
        m0=np.tile(m0_direction, (ntot, 1)),
        A0=A0,
        Ms=Ms,
        K0=K0,
        alpha=4000.0,
        gamma=0.0,
        cuda=cuda,
        cvode=cvode,
        useavgn=1,
    )
    problem.u_ea[:, :] = easy_axis[np.newaxis, :]
    problem.t = np.linspace(0.0, t_end, nt)
    problem.nt = len(problem.t)
    problem.t_conv = np.linspace(0.0, t_end, nt)
    problem.nt_conv = len(problem.t_conv)
    problem.exch_presize = problem.ntot * 28
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
    M_out = res[1][1, :, :, :]
    return tuple(
        Ms * np.mean(M_out[:, :n_steps, component], axis=0)
        for component in range(3)
    )


def output_stem(
    size_factor: float,
    n: int,
    use_fmm: bool,
    fmm_nterms: int,
    nlmax: int,
) -> str:
    """Return the legacy size-factor-based run identifier."""
    stem = f"single_grain_sf{size_factor:.2g}_n{n}"
    if use_fmm:
        stem += f"_fmm_on_N{fmm_nterms}_L{nlmax}"
    return stem


__all__ = [
    "DEFAULT_A0",
    "DEFAULT_K0",
    "DEFAULT_MS",
    "DEFAULT_MU0_MS_T",
    "DEFAULT_TILT_DEGREES",
    "MU0",
    "build_hysteresis_field",
    "characteristic_length",
    "create_single_grain_problem",
    "extract_mean_magnetisation",
    "interpolated_coercivity",
    "output_stem",
    "tilted_easy_axis",
    "tilted_field_direction",
    "z_easy_axis",
]
