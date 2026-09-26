"""Compare the field of a single inverted circular piece with a FEM simulation.

An inverted circular piece is the region between a circular arc and the straight
chord through its end points. Port of
matlab/examples/Magnetostatics/Validation_field_circpiece_inverted/MagTense_Validation_circpiece_inverted.m,
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


def validation_circpiece_inverted(show_plot: bool = True) -> list[float]:
    # center_pos is (r0, theta0, z0) and dev_center (dr, dtheta, dz)
    tile = Tiles(
        n=1,
        center_pos=[0.3, np.pi / 0.55, 0.6],
        dev_center=[0.15, np.pi / 6, 0.4],
        offset=[0.3, 0.5, 0.1],
        rot=[0, 0, 0],
        tile_type=4,
        M_rem=1.2 / MU0,
        easy_axis=[0.41562694, 0.41562694, 0.80901699],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )
    # Lines through the centre of the piece, evaluated at the FEM points
    files = [
        f"Validation_circpiece_inverted/Validation_circpiece_inverted_normH_{a}.txt"
        for a in "xyz"
    ]
    return norm_lines_validation(
        tile,
        "Inverted circular piece",
        files,
        [load_fem(f)[:, 0] for f in files],
        offset=[0.6271937452259475, 0.27251823835641853, 0.7],
        show_plot=show_plot,
    )


if __name__ == "__main__":
    validation_circpiece_inverted()
