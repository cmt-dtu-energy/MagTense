"""Smoke test of the dip-fmm demagnetization backend.

Runs a reduced version of fmm_vs_regular.py (one grid, a few order/depth pairs,
2 ns instead of 40 ns) and compares the final state of a relaxation with dip-fmm
against the dense demagnetization tensor, on the CPU and on the GPU. The error
has to fall with the expansion order. It takes about a minute on a 12^3 grid,
against the hours of the full sweep, so it is the check to run after a build.

    python fmm_smoke_test.py [cpu|cuda|both] [grid size, default 20]

The Python extension has to be built with USE_CDFMM=1.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fmm_vs_regular as example  # noqa: E402

# Keep the run short: relax for 2 ns instead of 40 ns.
example.SIMULATION_TIME = 2.0e-9
example.OUTPUT_STEPS = 11

which = sys.argv[1] if len(sys.argv) > 1 else "both"
modes = {"cpu": [False], "cuda": [True], "both": [False, True]}[which]
grid_size = int(sys.argv[2]) if len(sys.argv) > 2 else 20

failed = False
for use_cuda in modes:
    label = "CUDA" if use_cuda else "CPU"
    start = time.perf_counter()
    rows = example.run_sweep(
        grid_sizes=(grid_size,), orders=(2, 6, 10), depths=(2, 3), use_cuda=use_cuda
    )
    elapsed = time.perf_counter() - start
    print(f"\n{label}: order  depth  relative RMS   ({elapsed:.1f} s)")
    for row in rows:
        print(f"      {row['order']:>5}  {row['depth']:>5}  {row['relative_rms']:.3e}")
    for depth in (2, 3):
        errors = [r["relative_rms"] for r in rows if r["depth"] == depth]
        # Order 2 is a crude expansion; order 10 has to be far better and small.
        if not (errors[-1] < 1.0e-3 and errors[-1] < 0.1 * errors[0]):
            failed = True
            print(f"{label} depth {depth}: error does not converge: {errors}")

print("SMOKE_FAILED" if failed else "SMOKE_OK")
sys.exit(1 if failed else 0)
