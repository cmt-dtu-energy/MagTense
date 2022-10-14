#%%
import math
import matplotlib.pyplot as plt

from pathlib import Path

from magtense.magstatics import Tiles, run_simulation
from magtense.utils.eval import get_norm_magnetic_flux, load_COMSOL_eval, add_subplot


def valid_prism():
    eval_offset = [0.5, 0.4, 0.1]

    tile_val = Tiles(
        n=1,
        size=[0.6, 0.1, 0.3],
        offset=eval_offset,
        rot=[math.pi/2, -math.pi/3, math.pi/4],
        tile_type=2,
        M_rem=1.2/(4*math.pi*1e-7),
        easy_axis=[0.35355339, 0.61237244, 0.70710678],
        color=[1, 0, 0]
    )

    # Load reference points from COMSOL calculation
    COMSOL_eval_path = Path(__file__).parent.absolute() / '..' / '..' / '..' / \
        'documentation' / 'examples_FEM_validation' / 'Validation_prism'

    (pts_x, H_n_x_COMSOL) = load_COMSOL_eval('Validation_prism_normH_x.txt', eval_offset, COMSOL_eval_path)
    (pts_y, H_n_y_COMSOL) = load_COMSOL_eval('Validation_prism_normH_y.txt', eval_offset, COMSOL_eval_path)
    (pts_z, H_n_z_COMSOL) = load_COMSOL_eval('Validation_prism_normH_z.txt', eval_offset, COMSOL_eval_path)
    
    # x-axis
    (_, H_x) = run_simulation(tile_val, pts_x)
    H_n_x_MagTense  = get_norm_magnetic_flux(H_x)
    
    # y-axis
    (_, H_y) = run_simulation(tile_val, pts_y)
    H_n_y_MagTense  = get_norm_magnetic_flux(H_y)

    # z-axis
    (_, H_z) = run_simulation(tile_val, pts_z)
    H_n_z_MagTense  = get_norm_magnetic_flux(H_z)

    fig, ax = plt.subplots(1,3)
    fig.suptitle("PRISM - MagTensePython vs. COMSOL")
    add_subplot(ax[0], pts_x[:,0], 'x_axis', H_n_x_MagTense, H_n_x_COMSOL)
    add_subplot(ax[1], pts_y[:,1], 'y_axis', H_n_y_MagTense, H_n_y_COMSOL)
    add_subplot(ax[2], pts_z[:,2], 'z_axis', H_n_z_MagTense, H_n_z_COMSOL)
    plt.show()


if __name__ == '__main__':
    valid_prism()
# %%
