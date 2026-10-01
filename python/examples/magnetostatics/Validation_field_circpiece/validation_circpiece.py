"""Compare the field of a single permanent magnet circular piece with a FEM simulation.

Port of
matlab/examples/Magnetostatics/Validation_field_circpiece/MagTense_Validation_circpiece.m,
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


def validation_circpiece(show_plot: bool = True) -> list[float]:
    # center_pos is (r0, theta0, z0) and dev_center (dr, dtheta, dz): the piece spans
    # r0 +/- dr/2, theta0 +/- dtheta/2 and z0 +/- dz/2 around the tile's offset
    tile = Tiles(
        n=1,
        center_pos=[0.45, np.pi / 3, 0.35],
        dev_center=[0.5, np.pi / 7, 0.15],
        offset=[0.1, 0.3, 0.2],
        rot=[0, 0, 0],
        tile_type=3,
        M_rem=1.2 / MU0,
        easy_axis=[-0.3095974, -0.22493568, 0.92387953],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )
    # Lines through the centre of the piece, evaluated at the FEM points
    files = [f"Validation_circpiece/Validation_circpiece_normH_{a}.txt" for a in "xyz"]
    return norm_lines_validation(
        tile,
        "Circular piece",
        files,
        [load_fem(f)[:, 0] for f in files],
        offset=[0.406328622087633, 0.8631363102808783, 0.55],
        show_plot=show_plot,
    )


if __name__ == "__main__":
    validation_circpiece()
