"""Compare the field of a cone with a FEM simulation.

A cone tile is a right circular cone that stands on its base, and MagTense evaluates it
with the closed-form demagnetization tensor of R. Bjoerk, "The magnetic field of a
homogeneously magnetized cone". The cone has slant length 1.3 m and half-opening angle
0.6 rad, i.e. base radius 0.734 m and height 1.073 m, and carries the magnetization
(0.32, 0.74, 0.89) A/m. In the FEM model its apex is at the origin and its base in the
plane z = height, so the tile is turned upside down (rotation by pi about the x-axis) and
the centre of its base is moved to (0, 0, height). The field is compared with a COMSOL
solution along the line from (-0.6, -0.4, -0.1) to (1.1, 0.7, 1.3), which runs through
the magnet. The returned relative integrated errors in percent are those of H_x, H_y and
H_z along it.
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


def validation_cone(show_plot: bool = True) -> list[float]:
    slant, half_angle = 1.3, 0.6
    radius, height = slant * np.sin(half_angle), slant * np.cos(half_angle)
    M = np.array([0.32, 0.74, 0.89])
    line_start = np.array([-0.6, -0.4, -0.1])

    # size is (base radius, 0, height); the second entry is reserved for the top radius of
    # a truncated cone and must be zero. The tile has the centre of its base at offset and
    # its apex at the distance height along the rotated z-axis: turned upside down and
    # moved up by its height, the apex is at the origin and the base in the plane
    # z = height. With mu_r = 1 the magnetization is the remanence along the easy axis,
    # i.e. M.
    tile = Tiles(
        n=1,
        size=[radius, 0, height],
        offset=[0, 0, height],
        rot=[np.pi, 0, 0],
        tile_type=9,
        M_rem=np.linalg.norm(M),
        easy_axis=M / np.linalg.norm(M),
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )

    # Each FEM file holds one component of H in A/m against the x, y, z coordinates of the
    # mesh nodes on the line
    fem = [load_fem(f"Validation_cone/Validation_cone_{c}.txt") for c in ("Hx", "Hy", "Hz")]

    # In the FEM model the air around the cone is a sphere of radius 1.5 m about the middle
    # of the cone axis, enclosed by an infinite element domain. The line leaves the sphere
    # 2.4506 m from its start, 11 mm before its end. The coordinates are scaled in the
    # infinite element domain, so the nodes there cannot be compared and the comparison
    # stops at 2.45 m.
    line_length = 2.45

    errors, curves = [], []
    for i, comp in enumerate(("Hx", "Hy", "Hz")):
        # The export lists the nodes in mesh order, so they are sorted along the line.
        dist = np.linalg.norm(fem[i][:, :3] - line_start, axis=1)
        order = np.argsort(dist, kind="stable")
        order = order[dist[order] < line_length]
        fem_i, dist = fem[i][order], dist[order]
        pts = fem_i[:, :3]
        H = field(tile, pts)

        # A node on the surface of the cone is listed twice, with the field on each side
        # of the surface. The field jumps there, and the surface of the FEM mesh lies about
        # 0.2 mm inside the cone, so MagTense returns the field inside the cone for both.
        # These two nodes are left out of the error.
        twice = np.diff(dist) < 1e-12
        use = ~(np.append(twice, False) | np.append(False, twice))
        errors.append(
            calculate_relative_integral_error(dist[use], fem_i[use, 3], dist[use], H[use, i])
        )
        curves += [(f"MagTense, {comp}", dist, H[:, i], "rgb"[i] + "."),
                   (f"FEM, {comp}", dist, fem_i[:, 3], "rgb"[i] + "o")]

    report("Cone", errors, labels=("Hx", "Hy", "Hz"))
    if show_plot:
        plot_comparison(
            "Cone",
            [{"xlabel": "distance along the line from (-0.6, -0.4, -0.1) to (1.1, 0.7, 1.3) [m]",
              "curves": curves}],
            "H [A/m]",
        )
    return errors


if __name__ == "__main__":
    validation_cone()
