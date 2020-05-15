import os
import numpy as np
import matplotlib.pyplot as plt

import MagTense
from util_eval import load_COMSOL_eval, add_subplot

def main():
    # Define spherical tiles
    tile_val = MagTense.Tiles(1)

    # Define tile as sphere
    tile_val.set_tile_type(6)
    tile_val.set_size([0.25, 0, 0])
    tile_val.set_offset([0.3, 0.5, 0.1])
    tile_val.set_remanence(1.2 / (4*np.pi*1e-7))
    tile_val.set_easy_axis([[1, 0, -1]])

    eval_offset = [0.3, 0.4, 0.15]
    model_offset = [0.25, 0.25, 0.75]

    # Load reference points from COMSOL calculation
    COMSOL_eval_path = os.path.dirname(os.path.abspath(__file__)) \
        + '/../../documentation/examples_FEM_validation/Validation_sphere/'

    (eval_points_x, H_norm_x_COMSOL) = load_COMSOL_eval('py_Validation_sphere_normH_x.txt', \
        eval_offset, COMSOL_eval_path, model_offset=model_offset, unit='T')
    (eval_points_y, H_norm_y_COMSOL) = load_COMSOL_eval('py_Validation_sphere_normH_y.txt', \
        eval_offset, COMSOL_eval_path, model_offset=model_offset, unit='T')
    (eval_points_z, H_norm_z_COMSOL) = load_COMSOL_eval('py_Validation_sphere_normH_z.txt', \
        eval_offset, COMSOL_eval_path, model_offset=model_offset, unit='T')

    # x-axis
    (updated_tiles_x, H_x) = MagTense.run_simulation(tile_val, eval_points_x)
    H_norm_x_MagTense  = MagTense.get_norm_magnetic_flux(H_x)
    
    # y-axis
    (updated_tiles_y, H_y) = MagTense.run_simulation(tile_val, eval_points_y)
    H_norm_y_MagTense  = MagTense.get_norm_magnetic_flux(H_y)

    # z-axis
    (updated_tiles_z, H_z) = MagTense.run_simulation(tile_val, eval_points_z)
    H_norm_z_MagTense  = MagTense.get_norm_magnetic_flux(H_z)
    
    fig, ax = plt.subplots(1,3)
    fig.suptitle("SPHERE - MagTensePython vs. COMSOL")
    add_subplot(ax[0], eval_points_x[:,0], 'x_axis', H_norm_x_MagTense, H_norm_x_COMSOL)
    add_subplot(ax[1], eval_points_y[:,1], 'y_axis', H_norm_y_MagTense, H_norm_y_COMSOL)
    add_subplot(ax[2], eval_points_z[:,2], 'z_axis', H_norm_z_MagTense, H_norm_z_COMSOL)
    plt.show()

if __name__ == '__main__':
    main()
