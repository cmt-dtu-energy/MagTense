"""The field of a flat spiral coil, entered as a tile of type 101 (planar coil).

A planar coil tile is a flat coil in the xy-plane, centred on the tile's offset, made of
100 concentric circular loops evenly spaced from the inner radius a = size[0] to the
outer radius b = size[1]. Its M vector holds the current in each loop, in ampere, and
has to be the same in all three entries: M = [I, I, I]. The coil is not rotated and
does not take part in the magnetization iteration.

The field is checked against the Biot-Savart law summed over the same 100 loops, on the
axis, where each loop contributes I R^2 / (2 (R^2 + z^2)^(3/2)), and along a line off the
axis, where the loops are integrated numerically. Returns the maximum relative error
along each of the two lines.

Python counterpart of
matlab/examples/Magnetostatics/Example_006_planar_coil/MagTense_Example006_planar_coil.m.
"""

import numpy as np

from magtense.magstatics import Tiles, run_simulation

N_LOOPS = 100  # fixed in the Fortran implementation of the planar coil


def biot_savart_loops(radii: np.ndarray, current: float, pts: np.ndarray) -> np.ndarray:
    """H of coaxial circular loops in the xy-plane through the origin, by direct integration."""
    phi = np.linspace(0, 2 * np.pi, 4001)[:-1]
    dphi = phi[1] - phi[0]
    H = np.zeros_like(pts)
    for R in radii:
        # dl x (r - r') / |r - r'|^3 for every segment of the loop and every point
        src = np.c_[R * np.cos(phi), R * np.sin(phi), np.zeros_like(phi)]
        dl = np.c_[-R * np.sin(phi), R * np.cos(phi), np.zeros_like(phi)] * dphi
        d = pts[:, None, :] - src[None, :, :]
        dist3 = np.linalg.norm(d, axis=2) ** 3
        H += np.sum(np.cross(dl[None, :, :], d) / dist3[:, :, None], axis=1)
    return current * H / (4 * np.pi)


def planar_coil(
    a: float = 0.02, b: float = 0.05, current: float = 1.0, show_plot: bool = True
) -> list[float]:
    tiles = Tiles(n=1, size=[a, b, 0], offset=[0, 0, 0], tile_type=101, color=[0.8, 0.5, 0.2])
    tiles.M = ([current, current, current], 0)
    tiles.incl_it = (0, 0)

    radii = a + (b - a) * np.arange(N_LOOPS) / (N_LOOPS - 1)

    # On the axis, from the coil plane to four outer radii away
    z = np.linspace(0.001, 4 * b, 100)
    pts_axis = np.asfortranarray(np.c_[np.zeros_like(z), np.zeros_like(z), z])
    _, H_axis = run_simulation(tiles, pts_axis)
    Hz_exact = np.sum(
        current * radii[None, :] ** 2 / (2 * (radii[None, :] ** 2 + z[:, None] ** 2) ** 1.5),
        axis=1,
    )
    err_axis = np.max(np.abs(H_axis[:, 2] - Hz_exact) / Hz_exact)

    # Off the axis, along x at a height of one inner radius above the coil. A point in
    # the xz-plane is used because the planar coil is axisymmetric.
    x = np.linspace(0.0, 2 * b, 41)[1:]
    pts_off = np.asfortranarray(np.c_[x, np.zeros_like(x), np.full_like(x, a)])
    _, H_off = run_simulation(tiles, pts_off)
    H_off_exact = biot_savart_loops(radii, current, pts_off)
    err_off = np.max(
        np.linalg.norm(H_off - H_off_exact, axis=1) / np.linalg.norm(H_off_exact, axis=1)
    )

    print(f"Planar coil, a = {a} m, b = {b} m, {N_LOOPS} loops of {current} A")
    print(f"   max. relative error on the axis             {err_axis:.2e}")
    print(f"   max. relative error along x at z = a        {err_off:.2e}")

    if show_plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(z, H_axis[:, 2], "r.", label="MagTense")
        ax[0].plot(z, Hz_exact, "k-", label="Biot-Savart")
        ax[0].set_xlabel("z [m]")
        ax[0].set_ylabel(r"$H_z$ on the axis [A/m]")
        ax[1].plot(x, H_off[:, 0], "r.", label=r"MagTense, $H_x$")
        ax[1].plot(x, H_off[:, 2], "b.", label=r"MagTense, $H_z$")
        ax[1].plot(x, H_off_exact[:, 0], "r-", label=r"Biot-Savart, $H_x$")
        ax[1].plot(x, H_off_exact[:, 2], "b-", label=r"Biot-Savart, $H_z$")
        ax[1].set_xlabel(f"x [m], at z = {a} m")
        ax[1].set_ylabel("H [A/m]")
        ax[1].axvspan(a, b, color="0.9", label="coil windings")
        for axis in ax:
            axis.legend()
            axis.grid(True)
        fig.suptitle("Planar coil")
        fig.tight_layout()
        plt.show()

    return [float(err_axis), float(err_off)]


if __name__ == "__main__":
    planar_coil()
