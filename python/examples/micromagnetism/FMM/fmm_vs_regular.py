"""Compare dip-fmm and the regular MagTense demagnetisation calculation.

Build MagTense with ``USE_CDFMM=1`` before running this file. The ordinary
MagTense ``cuda`` setting controls both demagnetisation paths. Set
``MAGTENSE_USE_CUDA=0`` for a CPU-only smoke test.
"""

from pathlib import Path
import os
import sys

import numpy as np


# Allow the example to run directly from a MagTense source checkout.
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PYTHON_SOURCE = REPOSITORY_ROOT / "python" / "src"
if str(PYTHON_SOURCE) not in sys.path:
    sys.path.insert(0, str(PYTHON_SOURCE))

from magtense.micromag import MicromagProblem  # noqa: E402


GRID = (20, 20, 20)
CELL_SIZE = 5.0e-9
SATURATION_MAGNETISATION = 8.0e5
SIMULATION_TIME = 40.0e-9
OUTPUT_STEPS = 201
USE_CUDA = os.environ.get("MAGTENSE_USE_CUDA", "1").lower() not in {
    "0", "false", "no", "off"
}


def make_problem(*, use_dip_fmm: bool) -> MicromagProblem:
    """Return the same fixed micromagnetic problem with one demag backend changed."""
    cell_count = int(np.prod(GRID))
    initial_magnetisation = np.zeros((cell_count, 3))
    initial_magnetisation[:, 0] = 1.0
    initial_magnetisation[:, 2] = np.linspace(-0.1, 0.1, cell_count)
    initial_magnetisation /= np.linalg.norm(
        initial_magnetisation, axis=1, keepdims=True
    )

    problem = MicromagProblem(
        res=list(GRID),
        grid_L=np.asarray(GRID) * CELL_SIZE,
        grid_type="uniform",
        solver="dynamic",
        m0=initial_magnetisation,
        A0=1.3e-11,
        Ms=SATURATION_MAGNETISATION,
        K0=0.0,
        alpha=4.42e3,
        gamma=0.0,  # Pure relaxation keeps this long demonstration inexpensive.
        tol=1.0e-3,
        cuda=USE_CUDA,
        cvode=False,
        usereturnhall=False,
        use_cdfmm=use_dip_fmm,
        cdfmm_order=6,
        cdfmm_depth=3,
        cdfmm_basis="spherical",
    )
    problem.window_enabled = 0
    problem.trace_enabled = 0
    return problem


def zero_external_field(times: np.ndarray) -> np.ndarray:
    """Return a zero applied field for every requested time."""
    return np.zeros((len(times), 3))


def run_problem(*, use_dip_fmm: bool) -> tuple[np.ndarray, np.ndarray]:
    """Run 40 ns and return the output times and magnetisation trajectory."""
    problem = make_problem(use_dip_fmm=use_dip_fmm)
    result = problem.run_simulation(
        t_end=SIMULATION_TIME,
        nt=OUTPUT_STEPS,
        fct_h_ext=zero_external_field,
        nt_h_ext=2,
    )

    # M_out axes are (time, cell, applied-field index, Cartesian component).
    times = np.asarray(result[0])
    magnetisation = np.asarray(result[1][:, :, 0, :])
    return times, magnetisation


if __name__ == "__main__":
    # Regular MagTense call: dense demag tensor followed by its ordinary CPU path.
    print("\n=== Regular MagTense demag ===", flush=True)
    regular_times, regular_magnetisation = run_problem(use_dip_fmm=False)

    # FMM call: the physical problem is identical; only this flag changes.
    print("\n=== dip-fmm demag ===", flush=True)
    fmm_times, fmm_magnetisation = run_problem(use_dip_fmm=True)

    difference = fmm_magnetisation[-1] - regular_magnetisation[-1]
    relative_rms = np.linalg.norm(difference) / np.linalg.norm(
        regular_magnetisation[-1]
    )

    print(f"Cells: {int(np.prod(GRID))}")
    print(f"Simulated time: {SIMULATION_TIME * 1.0e9:.1f} ns")
    print(f"Output states: {len(regular_times)}")
    print(f"MagTense cuda setting: {USE_CUDA}")
    print(f"Final-state relative RMS difference: {relative_rms:.3e}")
