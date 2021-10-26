import os
import math
import matplotlib.pyplot as plt

import magtense
from util_eval import load_COMSOL_eval, add_subplot


def valid_circpiece():
    # Define validation tile
    tile_val = magtense.Tiles(1)
    tile_val.set_center_pos([[0.45, math.pi/3, 0.35]])
    tile_val.set_dev_center([0.5, math.pi/7, 0.15])
    tile_val.set_offset_i([0.1, 0.3, 0.2],0)
    tile_val.set_rotation_i([0, 0, 0],0)
    tile_val.set_tile_type(3)
    tile_val.set_remanence(1.2 / (4*math.pi*1e-7))
    tile_val.set_easy_axis([[-0.3095974 , -0.22493568,  0.92387953]])
    tile_val.set_color([[1, 0, 0]])

    eval_offset = [0.406328622087633, 0.8631363102808783, 0.55] 

    # Load reference points from COMSOL calculation
    COMSOL_eval_path = os.path.dirname(os.path.abspath(__file__)) + '/../../documentation/examples_FEM_validation/Validation_circpiece/'

    (eval_points_x, H_norm_x_COMSOL) = load_COMSOL_eval('Validation_circpiece_normH_x.txt', eval_offset, COMSOL_eval_path)
    (eval_points_y, H_norm_y_COMSOL) = load_COMSOL_eval('Validation_circpiece_normH_y.txt', eval_offset, COMSOL_eval_path)
    (eval_points_z, H_norm_z_COMSOL) = load_COMSOL_eval('Validation_circpiece_normH_z.txt', eval_offset, COMSOL_eval_path)

    # x-axis
    (updated_tiles_x, H_x) = magtense.run_simulation(tile_val, eval_points_x)
    H_norm_x_MagTense  = magtense.get_norm_magnetic_flux(H_x)
    
    # y-axis
    (updated_tiles_y, H_y) = magtense.run_simulation(tile_val, eval_points_y)
    H_norm_y_MagTense  = magtense.get_norm_magnetic_flux(H_y)

    # z-axis
    (updated_tiles_z, H_z) = magtense.run_simulation(tile_val, eval_points_z)
    H_norm_z_MagTense  = magtense.get_norm_magnetic_flux(H_z)
    
    fig, ax = plt.subplots(1,3)
    fig.suptitle("CIRCPIECE - MagTensePython vs. COMSOL")
    add_subplot(ax[0], eval_points_x[:,0], 'x_axis', H_norm_x_MagTense, H_norm_x_COMSOL)
    add_subplot(ax[1], eval_points_y[:,1], 'y_axis', H_norm_y_MagTense, H_norm_y_COMSOL)
    add_subplot(ax[2], eval_points_z[:,2], 'z_axis', H_norm_z_MagTense, H_norm_z_COMSOL)
    plt.show()


if __name__ == '__main__':
    valid_circpiece()
