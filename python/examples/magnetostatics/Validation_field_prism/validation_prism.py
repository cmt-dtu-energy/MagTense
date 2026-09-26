"""Compare the field of a single permanent magnet prism with a FEM simulation.

Port of matlab/examples/Magnetostatics/Validation_field_prism/MagTense_Validation_prism.m,
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


def validation_prism(show_plot: bool = True) -> list[float]:
    # The rotation is applied first around z, then around y and finally around x
    tile = Tiles(
        n=1,
        size=[0.6, 0.1, 0.3],
        offset=[0.5, 0.4, 0.1],
        rot=[np.pi / 2, -np.pi / 3, np.pi / 4],
        tile_type=2,
        M_rem=1.2 / MU0,
        easy_axis=[0.35355339, 0.61237244, 0.70710678],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )
    # Lines through the centre of the prism
    x = np.linspace(-0.5, 1.5, 2001)
    z = np.linspace(-1.0, 1.0, 2001)
    return norm_lines_validation(
        tile,
        "Prism",
        [f"Validation_prism/Validation_prism_normH_{a}.txt" for a in "xyz"],
        [x, x, z],
        offset=[0.5, 0.4, 0.1],
        show_plot=show_plot,
    )


if __name__ == "__main__":
    validation_prism()
