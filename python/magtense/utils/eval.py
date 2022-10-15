import numpy as np
import matplotlib.pyplot as plt

from typing import List, Optional
from pathlib import Path

from magtense.magstatics import Tiles, run_simulation


def get_average_magnetic_flux(H):
    mu0 = 4 * np.pi * 1e-7
    norm = [np.linalg.norm(H_point) * mu0 for H_point in H]
    return sum(norm) / len(norm)


def get_p2p(H):
    mu0 = 4 * np.pi * 1e-7
    norm = [np.linalg.norm(H_point) * mu0 for H_point in H]
    return max(norm) - min(norm)


def load_COMSOL(
    fname: str,
    eval_offset: List,
    COMSOL_eval_path: Path, 
    model_offset: List,
    unit: str,
    pts_special: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
    '''
    Load reference points from COMSOL calculation
    '''
    with open(Path(COMSOL_eval_path, fname), "r") as file:
        T = file.readlines()[8:]

    T_split = np.asarray([line.split() for line in T], dtype=np.float64)
    H_norm_COMSOL = T_split[:,1]
    if unit == 'T': H_norm_COMSOL *= 4 * np.pi * 1e-7
    pts_coor = T_split[:,0] if pts_special is None else pts_special
    struc = np.ones(len(pts_coor))

    if fname[-5] == 'x':
        pts = np.c_[pts_coor - model_offset[0],
                    struc * eval_offset[1],
                    struc * eval_offset[2]]
    elif fname[-5] == 'y':
        pts = np.c_[struc * eval_offset[0],
                    pts_coor - model_offset[1],
                    struc * eval_offset[2]]
    elif fname[-5] == 'z':
        pts = np.c_[struc * eval_offset[0],
                    struc * eval_offset[1],
                    pts_coor - model_offset[2]]

    return pts, H_norm_COMSOL


def validation(
    shape: str,
    tile: Tiles,
    offset: List,
    model_offset: List = [0, 0, 0],
    plot_COMSOL: bool = True,
    plot_error: bool = False,
    unit: str = 'A/m'
    ) -> None:
    mu0 = 4 * np.pi * 1e-7
    prefix = 'py_' if 'spher' in shape else ''
    suffix = '_prolate' if shape == 'spheroid' else ''
    COMSOL_eval_path = Path(__file__).parent.absolute() / '..' / '..' / '..' / \
        'documentation' / 'examples_FEM_validation' / f'Validation_{shape}'

    fig, ax = plt.subplots(1,3)
    fig.suptitle(f'{shape} - MagTensePy vs. COMSOL')

    for i, coord in enumerate(['x', 'y', 'z']):
        fname = f'{prefix}Validation_{shape}{suffix}_normH_{coord}.txt'
        pts, H_n_COMSOL = load_COMSOL(fname, offset, COMSOL_eval_path, model_offset, unit)
        _, H_mt = run_simulation(tile, pts)
        H_n_mt = [np.linalg.norm(H_point) * mu0 for H_point in H_mt]

        ax[i].plot(pts[:,i], H_n_mt, 'r*', label='MagTense')
        if plot_COMSOL: ax[i].plot(pts[:,i], H_n_COMSOL, 'bx', label='COMSOL')

        if plot_error:
            error = abs(H_n_COMSOL - H_n_mt)
            ax[i].plot(pts[:,i], error, 'g', label= 'Error')

        ax[i].legend()
        ax[i].set_xlabel(f'{coord}_axis')
        ax[i].set_ylabel('H_norm')

    plt.show()
