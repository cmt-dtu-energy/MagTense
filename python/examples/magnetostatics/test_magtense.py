from pathlib import Path

import numpy as np

from magtense.magstatics import Tiles, run_simulation
from magtense.utils import create_plot


def load_COMSOL(
    fname: str,
    eval_offset: list,
    COMSOL_eval_path: Path,
    model_offset: list,
    unit: str,
    pts_special: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load reference points from COMSOL calculation
    """
    with Path.open(Path(COMSOL_eval_path, fname), "r") as file:
        T = file.readlines()[8:]

    T_split = np.asarray([line.split() for line in T], dtype=np.float64)
    H_norm_COMSOL = T_split[:, 1]
    if unit == "T":
        H_norm_COMSOL *= 4 * np.pi * 1e-7
    pts_coor = T_split[:, 0] if pts_special is None else pts_special
    struc = np.ones(len(pts_coor))

    if fname[-5] == "x":
        pts = np.c_[
            pts_coor - model_offset[0], struc * eval_offset[1], struc * eval_offset[2]
        ]
    elif fname[-5] == "y":
        pts = np.c_[
            struc * eval_offset[0], pts_coor - model_offset[1], struc * eval_offset[2]
        ]
    elif fname[-5] == "z":
        pts = np.c_[
            struc * eval_offset[0], struc * eval_offset[1], pts_coor - model_offset[2]
        ]

    return pts, H_norm_COMSOL


def test_prism(
    shape: str = "prism", model_offset: tuple = (0, 0, 0), unit: str = ("A/m",)
) -> None:
    mu0 = 4 * np.pi * 1e-7
    tile = Tiles(
        n=1,
        size=[0.6, 0.1, 0.3],
        offset=[0.5, 0.4, 0.1],
        rot=[np.pi / 2, -np.pi / 3, np.pi / 4],
        tile_type=2,
        M_rem=1.2 / mu0,
        easy_axis=[0.35355339, 0.61237244, 0.70710678],
        color=[1, 0, 0],
    )
    offset = [0.5, 0.4, 0.1]

    mu0 = 4 * np.pi * 1e-7
    prefix = "py_" if "spher" in shape else ""
    suffix = "_prolate" if shape == "spheroid" else ""
    COMSOL_eval_path = (
        Path(__file__).parent.absolute()
        / ".."
        / ".."
        / ".."
        / "documentation"
        / "examples_FEM_validation"
        / f"Validation_{shape}"
    )

    for coord in ["x", "y", "z"]:
        fname = f"{prefix}Validation_{shape}{suffix}_normH_{coord}.txt"
        pts, H_n_COMSOL = load_COMSOL(
            fname, offset, COMSOL_eval_path, model_offset, unit
        )
        _, H_mt = run_simulation(tile, pts)
        H_n_mt = [np.linalg.norm(H_point) * mu0 for H_point in H_mt]

        print(f"Ten largest errors ({coord}): ", np.sort(abs(H_n_COMSOL - H_n_mt))[-5:])
        assert np.any(np.sort(abs(H_n_COMSOL - H_n_mt))[:-1] < 5e-3)


def test_plot_fn() -> None:
    """
    Test the plot function
    """
    mu0 = 4 * np.pi * 1e-7
    tiles = Tiles(
        n=6,
        M_rem=1.2 / mu0,
        tile_type=[2, 1, 3, 4, 5, 7],
        color=[
            [1, 0, 0],
            [0, 0, 1],
            [1, 0.5, 0],
            [0.3, 0.8, 0.2],
            [0, 0, 0],
            [1, 0, 1],
        ],
    )

    # 0: Prism
    tiles.size = ([0.1, 0.3, 0.2], 0)
    tiles.offset = ([0.1, 0.2, 0.1], 0)

    # 1: Cylindrical Tiles
    tiles.center_pos = ([1, 0, 0.3], 1)
    tiles.dev_center = ([0.15, np.pi / 9, 0.3], 1)

    # 2: Circpiece
    tiles.center_pos = ([0.85, np.pi / 5, 1.2], 2)
    tiles.dev_center = ([0.15, np.pi / 7, 0.25], 2)

    # 3: Inverted Circpiece
    tiles.center_pos = ([0.2, np.pi / 6, 0.75], 3)
    tiles.dev_center = ([0.05, np.pi / 4, 0.4], 3)

    # 4: Tetrahedron
    tiles.vertices = (
        np.array(
            [[0.65, 0.9, 0.5], [0.8, 0.9, 0.7], [0.85, 0.55, 0.25], [0.95, 0.85, 0.15]]
        ),
        4,
    )

    # 5: Prolate Spheroid
    tiles.size = ([0.1, 0.3, 0.1], 5)
    tiles.offset = ([0.1, 0.6, 0.7], 5)
    tiles.rot = ([0, 0, 2], 5)

    # Call the plot function
    create_plot(tiles, show=False)


def test_soft_sphere_in_uniform_field() -> None:
    """A soft sphere in a uniform applied field, entered as a tile of type 102.

    The sphere with constant permeability has the closed-form magnetization
    M = 3 (mu_r - 1) / (mu_r + 2) H_app, and the total field on the axis outside it is the
    applied field plus that of a point dipole with the sphere's moment.
    """
    mu0 = 4 * np.pi * 1e-7
    mu_r, radius = 20.0, 0.01
    H_app = np.array([0.0, 0.0, 0.1 / mu0])
    tiles = Tiles(
        n=1, size=[radius, 0, 0], offset=[0, 0, 0], tile_type=6, magnet_type=[3],
        mu_r_ea=mu_r, mu_r_oa=mu_r, M_rem=0.0, easy_axis=[0, 0, 1],
    )
    tiles.add_uniform_field(H_app)
    r = np.array([2.0, 5.0]) * radius
    pts = np.c_[np.zeros_like(r), np.zeros_like(r), r]
    tiles, H = run_simulation(tiles, pts, max_error=1e-10, max_it=500)

    M_exact = 3 * (mu_r - 1) / (mu_r + 2) * H_app
    V = 4 / 3 * np.pi * radius**3
    H_exact = H_app[None, :] + 2 * (M_exact * V)[None, :] / (4 * np.pi * r[:, None] ** 3)
    assert np.linalg.norm(tiles.M[0] - M_exact) / np.linalg.norm(M_exact) < 1e-6
    assert np.max(np.linalg.norm(H - H_exact, axis=1) / np.linalg.norm(H_exact, axis=1)) < 1e-6


def test_state_function_sphere_independent_of_mu_r() -> None:
    """A state-function tile is described by its curve alone.

    The internal field of a soft tile is solved with the tile's own material law, so the
    permeability entered for a state-function tile must not change the result, and the
    magnetization must satisfy M = M_curve(H_app - M/3) for a sphere. This guards against
    linearizing the self-consistency with mu_r_ea, which made the same sphere give anything
    from 0.06 T to 1.9 T depending on that parameter.
    """
    import importlib_resources

    from magtense.magstatics import iterate_magnetization

    mu0 = 4 * np.pi * 1e-7
    radius, H_app = 0.01, np.array([0.0, 0.0, 0.1 / mu0])
    M_norms = []
    for mu_r in (1.0, 5.0, 20.0, 100.0):
        tiles = Tiles(
            n=1, size=[radius, 0, 0], offset=[0, 0, 0], tile_type=6, magnet_type=[2],
            mu_r_ea=mu_r, mu_r_oa=mu_r, M_rem=0.0, easy_axis=[0, 0, 1],
        )
        tiles.stfcn_index = [1]
        tiles.add_uniform_field(H_app)
        tiles = iterate_magnetization(tiles, max_error=1e-10, max_it=500, mu_r=20)
        M_norms.append(np.linalg.norm(tiles.M[0]))
    M_norms = np.asarray(M_norms)
    assert np.max(np.abs(M_norms - M_norms[0])) / M_norms[0] < 1e-8

    # Self-consistency against the curve the solver used (linear interpolation of the table
    # against the solver's spline, hence the loose tolerance)
    ref = importlib_resources.files("magtense") / "mat/Fe_mur_20_Ms_2_1.csv"
    with importlib_resources.as_file(ref) as path:
        data = np.genfromtxt(path, delimiter=";", dtype=np.float64)
    H_int = np.linalg.norm(H_app) - M_norms[0] / 3
    M_curve = np.interp(H_int, data[1:, 0], data[1:, 1])
    assert abs(M_norms[0] - M_curve) / M_curve < 2e-2

