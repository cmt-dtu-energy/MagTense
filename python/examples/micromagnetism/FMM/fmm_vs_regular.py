"""Sweep dip-fmm order/depth and compare with regular MagTense demag.

The regular result is computed once per grid. Every dip-fmm run uses the same
physical problem and the persistent plans prepared by ``prepare_fmm_cache.py``.
Build MagTense with ``USE_CDFMM=1`` before running this file.
"""

import os
import sys
from pathlib import Path

import numpy as np


# Allow the example to run directly from a MagTense source checkout.
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PYTHON_SOURCE = REPOSITORY_ROOT / "python" / "src"
if str(PYTHON_SOURCE) not in sys.path:
    sys.path.insert(0, str(PYTHON_SOURCE))

from magtense.micromag import MicromagProblem  # noqa: E402


GRID_SIZES = (15, 20, 25, 30)
ORDERS = tuple(range(1, 11))
DEPTHS = (2, 3, 4, 5)

CELL_SIZE = 5.0e-9
SATURATION_MAGNETISATION = 8.0e5
SIMULATION_TIME = 40.0e-9
OUTPUT_STEPS = 201
USE_CUDA = os.environ.get("MAGTENSE_USE_CUDA", "1").lower() not in {
    "0", "false", "no", "off"
}


def make_problem(
    grid_size: int,
    *,
    use_dip_fmm: bool,
    order: int = 1,
    depth: int = 2,
    use_cuda: bool = USE_CUDA,
) -> MicromagProblem:
    """Create one cubic, non-periodic micromagnetic problem."""
    grid = (grid_size,) * 3
    cell_count = grid_size**3
    initial_magnetisation = np.zeros((cell_count, 3))
    initial_magnetisation[:, 0] = 1.0
    initial_magnetisation[:, 2] = np.linspace(-0.1, 0.1, cell_count)
    initial_magnetisation /= np.linalg.norm(
        initial_magnetisation, axis=1, keepdims=True
    )

    problem = MicromagProblem(
        res=list(grid),
        grid_L=np.asarray(grid) * CELL_SIZE,
        grid_type="uniform",
        solver="dynamic",
        m0=initial_magnetisation,
        A0=1.3e-11,
        Ms=SATURATION_MAGNETISATION,
        K0=0.0,
        alpha=4.42e3,
        gamma=0.0,  # Pure relaxation keeps the sweep inexpensive.
        tol=1.0e-3,
        cuda=use_cuda,
        cvode=False,
        usereturnhall=False,
        use_cdfmm=use_dip_fmm,
        cdfmm_order=order,
        cdfmm_depth=depth,
        cdfmm_basis="spherical",
    )
    problem.window_enabled = 0
    problem.trace_enabled = 0
    return problem


def zero_external_field(times: np.ndarray) -> np.ndarray:
    """Return a zero applied field for every requested time."""
    return np.zeros((len(times), 3))


def run_problem(
    grid_size: int,
    *,
    use_dip_fmm: bool,
    order: int = 1,
    depth: int = 2,
    use_cuda: bool = USE_CUDA,
) -> np.ndarray:
    """Run 40 ns and return only the final magnetisation state."""
    problem = make_problem(
        grid_size,
        use_dip_fmm=use_dip_fmm,
        order=order,
        depth=depth,
        use_cuda=use_cuda,
    )
    result = problem.run_simulation(
        t_end=SIMULATION_TIME,
        nt=OUTPUT_STEPS,
        fct_h_ext=zero_external_field,
        nt_h_ext=2,
    )
    # M_out axes are (time, cell, applied-field index, Cartesian component).
    return np.asarray(result[1][-1, :, 0, :]).copy()


def run_sweep(
    grid_sizes=GRID_SIZES,
    orders=ORDERS,
    depths=DEPTHS,
    *,
    use_cuda: bool = USE_CUDA,
) -> list[dict[str, int | float]]:
    """Run one regular reference per grid and every requested dip-fmm case."""
    rows = []
    for grid_size in grid_sizes:
        print(f"\n=== {grid_size}^3 regular reference ===", flush=True)
        regular_final = run_problem(
            grid_size, use_dip_fmm=False, use_cuda=use_cuda
        )

        for order in orders:
            for depth in depths:
                print(
                    f"\n=== {grid_size}^3, order={order}, depth={depth} ===",
                    flush=True,
                )
                fmm_final = run_problem(
                    grid_size,
                    use_dip_fmm=True,
                    order=order,
                    depth=depth,
                    use_cuda=use_cuda,
                )
                relative_rms = np.linalg.norm(
                    fmm_final - regular_final
                ) / np.linalg.norm(regular_final)
                row = {
                    "grid_size": grid_size,
                    "particles": grid_size**3,
                    "order": order,
                    "depth": depth,
                    "relative_rms": float(relative_rms),
                }
                rows.append(row)
                print(f"relative RMS: {relative_rms:.3e}", flush=True)
    return rows


if __name__ == "__main__":
    results = run_sweep()

    print("\nGrid     N  Order  Depth  Relative RMS")
    for row in results:
        print(
            f"{row['grid_size']:>2}^3  {row['particles']:>5}"
            f"  {row['order']:>5}  {row['depth']:>5}"
            f"  {row['relative_rms']:.3e}"
        )
