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
