"""Compare the field of a full cylinder with a FEM simulation.

A cylindrical tile that spans the full circle (dtheta = 2 pi) with inner radius zero is a
full cylinder, and MagTense evaluates it with the closed-form field of Caciagli et al.,
J. Magn. Magn. Mater. 456 (2018) 423, instead of the cylinder-piece integrals. The cylinder
has radius 1.1 m and height 0.75 m, is centred at the origin and carries the magnetization
(2, 3, 4) A/m, which has both an axial and a transverse component. The field is compared
with a COMSOL solution along the line from the origin to (-1.5, -0.75, 1.25), which runs
through the magnet and leaves it through the top end surface. The returned relative
integrated errors in percent are those of H_x, H_y and H_z along it.
"""

import sys
from pathlib import Path

import numpy as np

# fem_comparison lives one directory up, shared by all the validations
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fem_comparison import (  # noqa: E402
    calculate_relative_integral_error,
    field,
    load_fem,
    plot_comparison,
    report,
)

from magtense.magstatics import Tiles  # noqa: E402


def validation_cylinder(show_plot: bool = True) -> list[float]:
    radius, height = 1.1, 0.75
    M = np.array([2.0, 3.0, 4.0])

    # center_pos is (r0, theta0, z0) and dev_center (dr, dtheta, dz): the full circle with
    # inner radius r0 - dr/2 = 0 makes the tile a full cylinder. With mu_r = 1 the
    # magnetization is the remanence along the easy axis, i.e. M.
    tile = Tiles(
        n=1,
        center_pos=[radius / 2, 0, 0],
        dev_center=[radius, 2 * np.pi, height],
        offset=[0, 0, 0],
        rot=[0, 0, 0],
        tile_type=1,
        M_rem=np.linalg.norm(M),
        easy_axis=M / np.linalg.norm(M),
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )

    # Each FEM file holds one component of H in A/m against the x, y, z coordinates of the
    # evaluation points
    fem = [load_fem(f"Validation_cylinder/Validation_cylinder_full_{c}.txt")
           for c in ("Hx", "Hy", "Hz")]
    errors, curves = [], []
    for i, comp in enumerate(("Hx", "Hy", "Hz")):
        # The export lists the nodes in mesh order, and the node on the top end surface
        # twice, with the value from each side of the surface. The points are evaluated
        # as given (MagTense takes the limit from outside on a surface) and sorted along
        # the line for the integrated error.
        fem_i = fem[i][np.argsort(np.linalg.norm(fem[i][:, :3], axis=1), kind="stable")]
        pts = fem_i[:, :3]
        H = field(tile, pts)
        dist = np.linalg.norm(pts, axis=1)
        errors.append(calculate_relative_integral_error(dist, fem_i[:, 3], dist, H[:, i]))
        curves += [(f"MagTense, {comp}", dist, H[:, i], "rgb"[i] + "."),
                   (f"FEM, {comp}", dist, fem_i[:, 3], "rgb"[i] + "o")]

    report("Full cylinder", errors, labels=("Hx", "Hy", "Hz"))
    if show_plot:
        plot_comparison(
            "Full cylinder",
            [{"xlabel": "distance from the origin along the line to (-1.5, -0.75, 1.25) [m]",
              "curves": curves}],
            "H [A/m]",
        )
    return errors


if __name__ == "__main__":
    validation_cylinder()
