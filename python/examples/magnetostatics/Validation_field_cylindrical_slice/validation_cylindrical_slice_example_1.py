"""Compare the field of a cylindrical slice with a FEM simulation, example 1.

Port of
matlab/examples/Magnetostatics/Validation_field_cylindrical_slice/MagTense_Validation_cylindrical_slice_example_1.m,
with the same tile, evaluation points and error measure. The field is evaluated along a
straight line through space that starts at (2, -1, -3), and the returned relative
integrated errors in percent are those of H_x, H_y and H_z along it.
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


def validation_cylindrical_slice_example_1(show_plot: bool = True) -> list[float]:
    # center_pos is (r0, theta0, z0) and dev_center (dr, dtheta, dz)
    tile = Tiles(
        n=1,
        center_pos=[5.3984, np.pi / 8, 0],
        dev_center=[6.4672 - 4.3296, np.pi / 8, 1],
        offset=[0, 0, 0],
        rot=[0, 0, 0],
        tile_type=1,
        M_rem=1.2 / MU0,
        easy_axis=[1, 1, 1],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )

    # Each FEM file holds one component of H in A/m against one coordinate of the line
    fem = [load_fem(f"Validation_cylinder/Validation_cylinder_example_1_{c}.txt")
           for c in ("Hx", "Hy", "Hz")]
    # A small offset of 1e-3 keeps the points off the exact surface of the tile
    pts = np.column_stack([f[:, 0] + 1e-3 for f in fem])
    H = field(tile, pts)

    start = np.array([2, -1, -3])
    dist = np.linalg.norm(pts - start, axis=1)
    fem_dist = np.linalg.norm(np.column_stack([f[:, 0] for f in fem]) - start, axis=1)

    errors, curves = [], []
    for i, comp in enumerate(("Hx", "Hy", "Hz")):
        errors.append(calculate_relative_integral_error(fem_dist, fem[i][:, 1], dist, H[:, i]))
        curves += [(f"MagTense, {comp}", dist, H[:, i], "rgb"[i] + "."),
                   (f"FEM, {comp}", fem_dist, fem[i][:, 1], "rgb"[i] + "o")]

    report("Cylindrical slice, example 1", errors, labels=("Hx", "Hy", "Hz"))
    if show_plot:
        plot_comparison(
            "Cylindrical slice, example 1",
            [{"xlabel": "distance along the line [m]", "curves": curves}],
            "H [A/m]",
        )
    return errors


if __name__ == "__main__":
    validation_cylindrical_slice_example_1()
