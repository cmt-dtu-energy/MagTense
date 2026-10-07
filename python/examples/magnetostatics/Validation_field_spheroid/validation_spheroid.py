"""Compare the field of a single permanent magnet prolate spheroid with a FEM simulation.

Port of
matlab/examples/Magnetostatics/Validation_field_spheroid/MagTense_Validation_spheroid.m
(its prolate case, the one the MATLAB test suite runs), with the same tile, evaluation
points and error measure. Returns the relative integrated error in percent of H_x, H_y
and H_z along a line in x through the centre.
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


def _rot_z(t: float) -> np.ndarray:
    return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])


def _rot_y(t: float) -> np.ndarray:
    return np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])


def validation_spheroid(show_plot: bool = True) -> list[float]:
    # A prolate spheroid, a > b = c, whose symmetry axis points along axis_vector
    abc = [0.00611061, 0.00305531, 0.00305531]
    offset = [-0.05781690, 0.00804030, 0.01479210]
    axis_vector = np.array([0.101152, -0.0385283, -0.225353])

    # The direction of the symmetry axis, as getSpheroidRotAngles computes it in MATLAB
    # with rot_axis 'c': the local x-axis (the axis of revolution, as b = c) rotated by
    # the polar and azimuthal angles of axis_vector. For a spheroid the python Tiles take
    # this axis as rot and convert it to rotation angles in the same way as MATLAB.
    phi = np.arctan2(axis_vector[1], axis_vector[0])
    theta = np.arctan2(axis_vector[2], np.hypot(axis_vector[0], axis_vector[1]))
    symm_axis = _rot_z(phi) @ _rot_y(np.pi / 2 - theta) @ _rot_z(0.0) @ np.array([1.0, 0, 0])

    tile = Tiles(
        n=1,
        size=abc,
        offset=offset,
        rot=symm_axis,
        tile_type=7,
        M_rem=1.2 / MU0,
        easy_axis=[1.1, 0.5, 0.3],
        mu_r_ea=1.0,
        mu_r_oa=1.0,
        color=[1, 0, 0],
    )

    # A line along x through the centre of the spheroid; the FEM files hold H in A/m
    x = np.linspace(-0.13, 0.1, 2301)
    pts = np.column_stack([x, np.full_like(x, offset[1]), np.full_like(x, offset[2])])
    H = field(tile, pts)

    errors, curves = [], []
    for i, comp in enumerate(("Hx", "Hy", "Hz")):
        fem = load_fem(f"Validation_spheroid/Validation_spheroid_prolate_{comp}_x.txt")
        errors.append(calculate_relative_integral_error(fem[:, 0], fem[:, 1], x, H[:, i]))
        curves += [(f"MagTense, {comp}", x, H[:, i], "rgb"[i] + "."),
                   (f"FEM, {comp}", fem[:, 0], fem[:, 1], "rgb"[i] + "o")]

    report("Prolate spheroid", errors, labels=("Hx", "Hy", "Hz"))
    if show_plot:
        plot_comparison("Prolate spheroid", [{"xlabel": "x [m]", "curves": curves}], "H [A/m]")
    return errors


if __name__ == "__main__":
    validation_spheroid()
