"""A soft magnetic sphere in a uniform applied field.

The applied field is entered as a tile of type 102, which is not a geometry but a uniform
field source: its M vector holds the applied field H_app in A/m, it magnetizes the other tiles
in the iteration and it is included in the field returned at the evaluation points.

Two checks with closed-form answers:

1. A sphere with constant relative permeability mu_r magnetizes uniformly with
       M = 3 (mu_r - 1) / (mu_r + 2) H_app,
   because its demagnetization factor is 1/3 and the internal field is H_app - M/3.
2. Outside the sphere the field is the applied field plus that of a point dipole with the
   sphere's moment M V, so on the axis through the centre along H_app the returned total
   field is H_app + 2 M V / (4 pi r^3).

The same sphere with the state-function iron of MagTense is run as well, and its
magnetization is compared with the state function evaluated at the internal field
H_app - M/3, which is the self-consistency the solver imposes. The two agree up to the
difference between the solver's spline through the tabulated curve and the linear
interpolation used here; the permeability entered for the tile plays no role for a
state-function tile.
"""

import numpy as np

from magtense.magstatics import Tiles, iterate_magnetization, run_simulation


def soft_sphere_in_uniform_field(
    mu_r: float = 20.0,
    radius: float = 0.01,
    H_app_T: float = 0.1,
    verbose: bool = True,
) -> dict:
    """Run the sphere for a constant permeability and return the relative errors."""
    mu0 = 4 * np.pi * 1e-7
    H_app = np.array([0.0, 0.0, H_app_T / mu0])  # A/m along z

    tiles = Tiles(
        n=1,
        size=[radius, 0, 0],
        offset=[0, 0, 0],
        tile_type=6,            # sphere
        magnet_type=[3],        # soft, constant permeability
        mu_r_ea=mu_r,
        mu_r_oa=mu_r,
        M_rem=0.0,
        easy_axis=[0, 0, 1],
        color=[0, 0, 1],
    )
    tiles.add_uniform_field(H_app)

    # Total field on the axis, outside the sphere
    r = np.array([2, 3, 5, 10]) * radius
    pts = np.c_[np.zeros_like(r), np.zeros_like(r), r]

    tiles, H = run_simulation(tiles, pts, max_error=1e-10, max_it=500)
    M = tiles.M[0]

    M_exact = 3 * (mu_r - 1) / (mu_r + 2) * H_app
    V = 4 / 3 * np.pi * radius**3
    H_exact = H_app[None, :] + 2 * (M_exact * V)[None, :] / (4 * np.pi * r[:, None] ** 3)

    err_M = np.linalg.norm(M - M_exact) / np.linalg.norm(M_exact)
    err_H = np.max(np.linalg.norm(H - H_exact, axis=1) / np.linalg.norm(H_exact, axis=1))

    if verbose:
        print(f"mu0 H_app = {H_app_T} T along z, sphere of radius {radius} m, mu_r = {mu_r}")
        print(f"   M        = {mu0 * M} T, exact {mu0 * M_exact} T, relative error {err_M:.2e}")
        print(f"   total field on the axis at r/a = {r / radius}: max relative error {err_H:.2e}")

    return {"M": M, "M_exact": M_exact, "err_M": err_M, "err_H": err_H, "H": H, "H_exact": H_exact}


def soft_sphere_state_function(radius: float = 0.01, H_app_T: float = 0.1, verbose: bool = True) -> float:
    """The same sphere with the state-function iron; returns the self-consistency error."""
    mu0 = 4 * np.pi * 1e-7
    H_app = np.array([0.0, 0.0, H_app_T / mu0])

    tiles = Tiles(
        n=1,
        size=[radius, 0, 0],
        offset=[0, 0, 0],
        tile_type=6,
        magnet_type=[2],        # soft, state function
        mu_r_ea=20,
        mu_r_oa=20,
        M_rem=0.0,
        easy_axis=[0, 0, 1],
        color=[0, 0, 1],
    )
    tiles.stfcn_index = [1]
    tiles.add_uniform_field(H_app)
    tiles = iterate_magnetization(tiles, max_error=1e-10, max_it=500, mu_r=20)
    M = tiles.M[0]

    # The state function M(H) that iterate_magnetization used, from the package data
    import importlib_resources

    ref = importlib_resources.files("magtense") / "mat/Fe_mur_20_Ms_2_1.csv"
    with importlib_resources.as_file(ref) as path:
        data = np.genfromtxt(path, delimiter=";", dtype=np.float64)
    H_tab, M_tab = data[1:, 0], data[1:, 1]

    # Internal field of a uniformly magnetized sphere, and the magnetization it should give
    H_int = np.linalg.norm(H_app - M / 3)
    M_should = np.interp(H_int, H_tab, M_tab)
    err = abs(np.linalg.norm(M) - M_should) / M_should

    if verbose:
        print(f"State-function iron: mu0 M = {mu0 * np.linalg.norm(M):.4f} T, "
              f"mu0 H_int = {mu0 * H_int:.4f} T, M(H_int) from the table {mu0 * M_should:.4f} T, "
              f"relative error {err:.2e}")
    return err


if __name__ == "__main__":
    soft_sphere_in_uniform_field()
    soft_sphere_state_function()
