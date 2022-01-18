import os
import math
import matplotlib.pyplot as plt

from magtense import magtense
from magtense.utils.eval import get_norm_magnetic_flux, load_COMSOL_eval, add_subplot


def valid_circ_inv():    
    # Define validation tile
    tile_val = magtense.Tiles(1)
    tile_val.set_center_pos([[0.3, math.pi/0.55, 0.6]])
    tile_val.set_dev_center([0.15, math.pi/6, 0.4])
    tile_val.set_offset_i([0.3, 0.5, 0.1],0)
    tile_val.set_rotation_i([0, 0, 0],0)
    tile_val.set_tile_type(4)
    tile_val.set_remanence(1.2 / (4*math.pi*1e-7))
    tile_val.set_easy_axis([[0.41562694, 0.41562694, 0.80901699]])
    tile_val.set_color([[1, 0, 0]])

    eval_offset = [0.6271937452259475, 0.27251823835641853, 0.7]

    # Load reference points from COMSOL calculation
    COMSOL_eval_path = os.path.dirname(os.path.abspath(__file__)) \
        + '/../../documentation/examples_FEM_validation/Validation_circpiece_inverted/'

    (eval_points_x, H_norm_x_COMSOL) = load_COMSOL_eval(
        'Validation_circpiece_inverted_normH_x.txt', eval_offset, COMSOL_eval_path)
    (eval_points_y, H_norm_y_COMSOL) = load_COMSOL_eval(
        'Validation_circpiece_inverted_normH_y.txt', eval_offset, COMSOL_eval_path)
    (eval_points_z, H_norm_z_COMSOL) = load_COMSOL_eval(
        'Validation_circpiece_inverted_normH_z.txt', eval_offset, COMSOL_eval_path)

    # x-axis
    (updated_tiles_x, H_x) = magtense.run_simulation(tile_val, eval_points_x)
    H_norm_x_MagTense  = get_norm_magnetic_flux(H_x)
    
    # y-axis
    (updated_tiles_y, H_y) = magtense.run_simulation(tile_val, eval_points_y)
    H_norm_y_MagTense  = get_norm_magnetic_flux(H_y)

    # z-axis
    (updated_tiles_z, H_z) = magtense.run_simulation(tile_val, eval_points_z)
    H_norm_z_MagTense  = get_norm_magnetic_flux(H_z)
    
    fig, ax = plt.subplots(1,3)
    fig.suptitle("CIRCPIECE_INV - MagTensePython vs. COMSOL")
    add_subplot(ax[0], eval_points_x[:,0], 'x_axis', H_norm_x_MagTense, H_norm_x_COMSOL)
    add_subplot(ax[1], eval_points_y[:,1], 'y_axis', H_norm_y_MagTense, H_norm_y_COMSOL)
    add_subplot(ax[2], eval_points_z[:,2], 'z_axis', H_norm_z_MagTense, H_norm_z_COMSOL)
    plt.show()


if __name__ == '__main__':
    valid_circ_inv()
