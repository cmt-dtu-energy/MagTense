"""Compare the field of a single permanent magnet sphere with a FEM simulation.

Port of matlab/examples/Magnetostatics/Validation_field_sphere/MagTense_Validation_sphere.m,
with the same tile, evaluation points and error measure. Returns the relative
integrated error in percent of H_x, H_y and H_z along a line in x through the centre.
"""

import sys
from pathlib import Path

import numpy as np

# fem_comparison lives one directory up, shared by all the validations
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fem_comparison import (  # noqa: E402
    MU0,
    calculate_relative_integral_error,
    field,
    load_fem,
    plot_comparison,
    report,
)

from magtense.magstatics import Tiles  # noqa: E402


def validation_sphere(show_plot: bool = True) -> list[float]:
    offset = [2, 3, 4]
    # size[0] is the radius
    tile = Tiles(
        n=1,
        size=[1.6, 0, 0],
        offset=offset,
        tile_type=6,
        M_rem=1.2 / MU0,
        easy_axis=[1.1, 0.5, 0.3],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )

    # A line along x through the centre of the sphere; the FEM files hold H in A/m
    x = np.linspace(-10, 10, 2001)
    pts = np.column_stack([x, np.full_like(x, offset[1]), np.full_like(x, offset[2])])
    H = field(tile, pts)

    errors, curves = [], []
    for i, comp in enumerate(("Hx", "Hy", "Hz")):
        fem = load_fem(f"Validation_sphere/Validation_sphere_{comp}_x.txt")
        errors.append(calculate_relative_integral_error(fem[:, 0], fem[:, 1], x, H[:, i]))
        curves += [(f"MagTense, {comp}", x, H[:, i], "rgb"[i] + "."),
                   (f"FEM, {comp}", fem[:, 0], fem[:, 1], "rgb"[i] + "o")]

    report("Sphere", errors, labels=("Hx", "Hy", "Hz"))
    if show_plot:
        plot_comparison("Sphere", [{"xlabel": "x [m]", "curves": curves}], "H [A/m]")
    return errors


if __name__ == "__main__":
    validation_sphere()
