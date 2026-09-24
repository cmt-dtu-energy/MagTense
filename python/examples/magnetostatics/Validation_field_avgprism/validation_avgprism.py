"""Compare the volume-averaged field of a prism with a FEM simulation.

The averaged prism (tile type 8) returns, at each evaluation point, the field of the
prism averaged over a rectangular observation volume centred on that point, whose size
is passed as obs_size. This is what a finite-size sensor, or a micromagnetic cell,
sees. Python counterpart of
matlab/examples/Magnetostatics/Validation_field_avgprism/MagTense_Validation_avgprism.m.
Returns the relative integrated error in percent of the x and y components.
"""

import sys
from pathlib import Path

import numpy as np

# fem_comparison lives one directory up, shared by all the validations
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fem_comparison import MU0, calculate_relative_integral_error, load_fem  # noqa: E402

from magtense.magstatics import Tiles, run_simulation  # noqa: E402


def validation_avgprism(show_plot: bool = True) -> list[float]:
    # A 3 x 5 x 1 prism at the origin, magnetized with 1 A/m along x
    tile = Tiles(
        n=1,
        size=[3.0, 5.0, 1.0],
        offset=[0.0, 0.0, 0.0],
        tile_type=8,
        M_rem=1.0,
        easy_axis=[1, 0, 0],
        color=[1, 0, 0],
    )

    # The COMSOL reference: the x coordinate of the observation volume and the average
    # of B_x and B_y over it, in tesla
    data = load_fem("Validation_avgprism/Avg_validation_comsol.txt")
    x = data[:, 0]

    # Observation volumes of 2 x 4 x 3 centred on a line along x at y = z = 5
    pts = np.asfortranarray(np.column_stack([x, np.full_like(x, 5.0), np.full_like(x, 5.0)]))
    obs_size = np.asfortranarray(np.tile([2.0, 4.0, 3.0], (len(x), 1)))
    _, H = run_simulation(tile, pts, obs_size=obs_size)
    B = H * MU0

    errors = [calculate_relative_integral_error(x, data[:, i + 1], x, B[:, i]) for i in range(2)]
    print(
        "Averaged prism: relative integrated error between MagTense and FEM is "
        f"<B_x> = {errors[0]:.3g} %, <B_y> = {errors[1]:.3g} %"
    )

    if show_plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 5))
        plt.plot(x, data[:, 1], "ko", fillstyle="none", label=r"$\langle B_x \rangle$, FEM")
        plt.plot(x, data[:, 2], "ro", fillstyle="none", label=r"$\langle B_y \rangle$, FEM")
        plt.plot(x, B[:, 0], "k", label=r"$\langle B_x \rangle$, MagTense")
        plt.plot(x, B[:, 1], "r", label=r"$\langle B_y \rangle$, MagTense")
        plt.xlabel("x [m]")
        plt.ylabel("Field averaged over the observation volume [T]")
        plt.title("Averaged prism - MagTense vs. FEM")
        plt.legend()
        plt.grid(True)
        plt.show()

    return errors


if __name__ == "__main__":
    validation_avgprism()
