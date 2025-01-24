import numpy as np
from magtense.magstatics import Tiles, run_simulation

from pathlib import Path
from typing import Optional, List


def load_COMSOL(
    fname: str,
    eval_offset: List,
    COMSOL_eval_path: Path,
    model_offset: List,
    unit: str,
    pts_special: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load reference points from COMSOL calculation
    """
    with open(Path(COMSOL_eval_path, fname), "r") as file:
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


def test_prism():
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
    shape = "prism"
    offset = [0.5, 0.4, 0.1]
    model_offset = ([0, 0, 0],)
    unit: str = ("A/m",)

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

    for i, coord in enumerate(["x", "y", "z"]):
        fname = f"{prefix}Validation_{shape}{suffix}_normH_{coord}.txt"
        pts, H_n_COMSOL = load_COMSOL(
            fname, offset, COMSOL_eval_path, model_offset, unit
        )
        _, H_mt = run_simulation(tile, pts)
        H_n_mt = [np.linalg.norm(H_point) * mu0 for H_point in H_mt]
        assert abs(H_n_COMSOL - H_n_mt) < 1e-6
