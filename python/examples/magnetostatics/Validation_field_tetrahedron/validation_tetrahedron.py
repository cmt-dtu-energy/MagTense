"""Compare the field of a single permanent magnet tetrahedron with a FEM simulation.

Port of
matlab/examples/Magnetostatics/Validation_field_tetrahedron/MagTense_Validation_tetrahedron.m,
with the same tile, evaluation points and error measure. Returns the relative
integrated error in percent of |mu0 H| along x, y and z.
"""

import sys
from pathlib import Path

import numpy as np

# fem_comparison lives one directory up, shared by all the validations
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fem_comparison import MU0, norm_lines_validation  # noqa: E402

from magtense.magstatics import Tiles  # noqa: E402


def validation_tetrahedron(show_plot: bool = True) -> list[float]:
    # A tetrahedron is given by its four vertices; offset and rot are not used.
    # |M| is 1.2 T / mu0 along the easy axis.
    tile = Tiles(
        n=1,
        vertices=[[2.5, 3, 1], [2, 1, 4], [1.5, 4, 3], [4.5, 5, 2]],
        tile_type=5,
        M_rem=1.2 / MU0,
        easy_axis=[0.324264068, 0.734846928, 0.891545179],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )
    x = np.linspace(-10, 10, 2001)
    return norm_lines_validation(
        tile,
        "Tetrahedron",
        [f"Validation_tetrahedron/Validation_tetrahedron_normH_{a}.txt" for a in "xyz"],
        [x, x, x],
        offset=[3, 3, 2.5],
        show_plot=show_plot,
    )


if __name__ == "__main__":
    validation_tetrahedron()
