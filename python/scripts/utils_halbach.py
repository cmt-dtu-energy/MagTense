#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from magtense.halbach import HalbachCylinder, EvaluationPoints
from magtense.magstatics import Tiles, iterate_magnetization, get_demag_tensor, get_H_field
from magtense.utils import plot_cube
from typing import Optional, List, Union
from pathlib import Path


def load_demag_tensor(
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
    store: bool = False
) -> np.ndarray:
    '''
    If store is set, then pre-calculation of the demagnetization tensor,
    which speeds up subsequent field calculation.
    Then, H = N * M becomes only a matrix multiplication.
    '''
    fname = f'N_{pts_eval.res[0]}_{pts_eval.res[1]}_{pts_eval.res[2]}_r_{pts_eval.r}_' \
            + f'hard_{halbach.n_layers}_{halbach.n_segs}_{halbach.n_axial}_' \
            + f'shim_{halbach.shim_layers}_{halbach.shim_segs}.npy'

    fpath = Path(__file__).parent.absolute() / '..' / 'demag' / fname

    if Path(fpath).is_file():
        demag_tensor = np.load(fpath)
    else:
        if not store:
            print('[WARNING] No demagnetization tensor found for this Halbach configuration!')
            print('[RECOMMENDATION] Run `load_demag_tensor(store=True)` first for a larger amount of samples.')

        print('Calculating demagnetization tensor...')
        demag_tensor = get_demag_tensor(halbach.tiles, pts_eval.coords)

        if store:
            if not fpath.parent.exists(): fpath.parent.mkdir()
            np.save(fpath, demag_tensor)

    return demag_tensor


def calc_demag_field(
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
) -> np.ndarray:
    '''
    Calculation of the resulting demagnetizing field inside the Halbach cylinder bore.
    MagTense uses herefore three steps:
    (1) Demagnetization tensor N
        This might be pre-calculated from the existig structure and a full shim magnet matrix.
    (2) Iterating the magnetizations M of each single magnetic tile.
        This becomes especially computationally expensive, when a lof of shim magnets are present.
    (3) Demagnetizing field: H = N * M
    '''
    demag_tensor = load_demag_tensor(halbach, pts_eval)
    it_tiles = iterate_magnetization(halbach.tiles, max_error=1e-12, mu_r=halbach.mu_r)
    field = get_H_field(it_tiles, pts_eval.coords, demag_tensor)
    # Filter posssible nan values - TODO Check with magtense
    idx_nan = [i for i, H_vec in enumerate(field) if np.isnan(np.linalg.norm(H_vec))]
    field[idx_nan,:] = [0,0,0]

    return field


def get_shim_mat(
    mat: Union[List, np.ndarray],
    layers: int,
    segs: int,
    symm: bool
) -> np.ndarray:
    '''
    Args:
        mat: Shim matrix solution.
        layers: Number of shim layers.
        segs: Shim magnets per layer.
        symm: If set, symmetry along center plane is assumed.

    Returns:
        Full shim matrix.
    '''
    n_shim = layers * segs
    shim_mat = np.zeros(n_shim)
    if symm:
        # Set bottom and center layers
        shim_mat[:segs * (layers // 2 + layers % 2)] = mat

        # Set top layers symmetrical to bottom layers
        for i in range(layers // 2):
            shim_mat[n_shim-(i+1)*segs:n_shim-i*segs] = mat[i*segs:(i+1)*segs]
    else:
        shim_mat[:] = mat

    return shim_mat


def eval_shimming(
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
    log_p2p: bool = False
) -> tuple[np.ndarray, float]:
    '''
    Args:
        halbach: Halbach cylinder, which produces magnetic field.
        pts_eval: Coordinates of evaluation points.
        log_p2p: If set, then log10(p2p) is returned.

    Returns:
        (1) Magnetic flux density.
        (2) Peak-to-peak (p2p) value in [parts per million].
    '''
    mu0 = 4 * np.pi * 1e-7
    h = calc_demag_field(halbach, pts_eval)
    b = h.reshape((*pts_eval.res,3)).transpose((3,0,1,2)) * mu0

    b_norm = [np.linalg.norm(h_vec) * mu0
              for i, h_vec in enumerate(h) 
              if (pts_eval.coords[i][0]**2 \
                + pts_eval.coords[i][1]**2 \
                + pts_eval.coords[i][2]**2) <= pts_eval.r**2]
    
    p2p = (max(b_norm) - min(b_norm)) / max(b_norm)
    p2p_res = np.log10(p2p) if log_p2p else p2p * 1e6

    return b, p2p_res


### API ###
def run_halbach_environment(
    solution: np.ndarray,
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
    symm: bool,
    log_p2p: bool = False,
):
    shim_mat = get_shim_mat(solution, halbach.shim_layers, halbach.shim_segs, symm)
    halbach.set_shim_matrix(shim_mat)
    b, p2p = eval_shimming(halbach, pts_eval, log_p2p)
    
    return b, p2p


### PLOTTING ###
def plot_shim_matrix(
    shim_mat: Union[List[int], np.ndarray],
    pts: np.ndarray,
    n_layers: int,
    n_segs: int,
) -> None:
    n_fe = pts.shape[0]
    r = 0.02

    if len(shim_mat) != n_fe:
        mat = np.zeros(n_fe)
        mat[:n_segs * (n_layers // 2 + n_layers % 2)] = shim_mat
        for i in range(n_layers // 2):
            mat[n_fe-((i + 1) * n_segs):n_fe-(i * n_segs)] \
                = shim_mat[i * n_segs:(i + 1) * n_segs]
        shim_mat = mat

    shim_magnets = Tiles(n=np.count_nonzero(shim_mat))
    idx = 0
    for shim_idx in range(n_fe):
        z = shim_idx // n_segs
        if shim_mat[shim_idx] == 1:
            shim_magnets.offset = ([pts[shim_idx, 0] * (z/2 + 1),
                                    pts[shim_idx, 1] * (z/2 + 1),
                                    pts[shim_idx, 2]], idx)
            shim_magnets.tile_type = (2, idx)
            shim_magnets.size = ([r / n_segs, r / n_segs, r], idx)
            color = [0.4, 0.4, 0.4]
            color[z - 3 * (z // 3)] = 0.9
            shim_magnets.color = (color, idx)
            idx += 1
        
        shim_idx += 1

    ax = plt.figure().gca(projection='3d')
    for i in range(shim_magnets.n):
        plot_cube(ax, shim_magnets.size[i], shim_magnets.offset[i],
                  shim_magnets.rot[i], shim_magnets.M[i],
                  shim_magnets.color[i])
    
    ax.set_xlim([-r, r])
    ax.set_ylim([-r, r])
    ax.set_zlim([-2 * r, 2 * r])
    ax.view_init(elev=90, azim=-90)
    plt.axis('off')
    plt.show()


def plot_field(field: np.ndarray) -> None:  
    z_idx = 0 if field.shape[0] == 2 else field.shape[3] // 2
    field = np.linalg.norm(field, axis=0)[:,:,z_idx]
    var_field = (field - field.mean()) * 1e6
    norm = colors.Normalize(vmin=var_field.min(), vmax=var_field.max())
    plt.close('all')
    ax = plt.gca()
    im = ax.imshow(var_field, cmap='rainbow', norm=norm, origin="lower")
    ax.set_title(f'Variation of |B| [ppm] - Avg: {field.mean():.3f} T')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig = plt.gcf()
    cbar_ax = fig.add_axes([0.825, 0.345, 0.015, 0.3])
    fig.colorbar(im, cax=cbar_ax)
    plt.show()


def plot_field_line(
    field: np.ndarray,
    eval_pts: Union[List[float], np.ndarray],
    best_p2p: Optional[int] = None
) -> None:
    res = field.shape[1] // 3
    field = np.linalg.norm(field, axis=0)
    plt.close('all')
    plt.plot(eval_pts, field[:res], color='r')
    plt.plot(eval_pts, field[res:2*res], color='g')
    plt.plot(eval_pts, field[2*res:], color='b')
    plt.legend(['x', 'y', 'z'])
    if best_p2p: plt.title(f'p2p: {best_p2p:.0f} ppm')
    plt.show()
#%%
