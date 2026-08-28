"""Prepare the dip-fmm plans used by the MagTense comparison sweep.

Existing validated cache entries are loaded instead of being rebuilt. The
plans are non-periodic, use the spherical basis and exactly the same uniform
cuboid geometry constructor as a normal MagTense simulation.
"""

import os

from fmm_vs_regular import (
    DEPTHS,
    GRID_SIZES,
    ORDERS,
    make_problem,
    zero_external_field,
)


# Disk plans are backend-independent, so avoid needless GPU uploads by default.
# MAGTENSE_USE_CUDA=1 can be used to exercise the production CUDA setup path.
USE_CUDA = os.environ.get("MAGTENSE_USE_CUDA", "0").lower() not in {
    "0", "false", "no", "off"
}


def prepare_plan(grid_size: int, order: int, depth: int) -> None:
    """Initialize one plan and skip the micromagnetic time evolution."""
    problem = make_problem(
        grid_size,
        use_dip_fmm=True,
        order=order,
        depth=depth,
        use_cuda=USE_CUDA,
    )
    problem.dummy_run = 1
    problem.run_simulation(
        t_end=1.0e-12,
        nt=2,
        fct_h_ext=zero_external_field,
        nt_h_ext=2,
    )


if __name__ == "__main__":
    total = len(GRID_SIZES) * len(ORDERS) * len(DEPTHS)
    case = 0
    for grid_size in GRID_SIZES:
        for order in ORDERS:
            for depth in DEPTHS:
                case += 1
                print(
                    f"\n[{case}/{total}] {grid_size}^3, order={order}, depth={depth}",
                    flush=True,
                )
                prepare_plan(grid_size, order, depth)

    print(f"\nPrepared {total} non-periodic spherical plan configurations.")
