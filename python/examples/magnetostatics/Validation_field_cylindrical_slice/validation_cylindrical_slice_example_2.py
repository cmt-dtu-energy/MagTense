"""Compare the field of a cylindrical slice with a FEM simulation, example 2.

Port of
matlab/examples/Magnetostatics/Validation_field_cylindrical_slice/MagTense_Validation_cylindrical_slice_example_2.m,
with the same tile, evaluation points and error measure. Returns the relative
integrated error in percent of |mu0 H| along x, y and z.
"""

import sys
from pathlib import Path

import numpy as np

# fem_comparison lives one directory up, shared by all the validations
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fem_comparison import MU0, load_fem, norm_lines_validation  # noqa: E402

from magtense.magstatics import Tiles  # noqa: E402


def validation_cylindrical_slice_example_2(show_plot: bool = True) -> list[float]:
    # center_pos is (r0, theta0, z0) and dev_center (dr, dtheta, dz)
    tile = Tiles(
        n=1,
        center_pos=[0.3, np.pi / 2, 0.5],
        dev_center=[0.3, np.pi / 4, 0.1],
        offset=[0.8, -0.1, 0.3],
        rot=[0, 0, 0],
        tile_type=1,
        M_rem=1.2 / MU0,
        easy_axis=[0.35355339, 0.35355339, 0.8660254],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )
    # Lines through (0.8, 0.2, 0.8) at the FEM points, shifted by 1e-6 to keep them off
    # the exact faces of the tile
    files = [f"Validation_cylinder/Validation_cylinder_example_2_normH_{a}.txt" for a in "xyz"]
    return norm_lines_validation(
        tile,
        "Cylindrical slice, example 2",
        files,
        [load_fem(f)[:, 0] + 1e-6 for f in files],
        offset=[0.8, 0.2, 0.8],
        show_plot=show_plot,
    )


if __name__ == "__main__":
    validation_cylindrical_slice_example_2()
