import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../source')
import MagTense
import MagTenseStandalone
import util_plot

def main():
    # Define validation tile
    offset = [0.5, 0.4, 0.1]

    tile_val = MagTense.Tiles(1)
    tile_val.set_size([[0.6, 0.1, 0.3]])
    tile_val.set_offset_i(offset,0)
    tile_val.set_rotation_i([math.pi/2, -math.pi/3, math.pi/4],0)
    tile_val.set_tile_type(2)
    tile_val.set_remanence(1.2 / (4*math.pi*1e-7))
    tile_val.set_mag_angle([[math.pi/3, math.pi/4]])
    tile_val.set_color([[1, 0, 0]])

    tile_standalone = MagTenseStandalone.Tile()
    tile_standalone.set_size([0.6, 0.1, 0.3]) 
    tile_standalone.set_offset(offset)
    tile_standalone.set_rotation([math.pi/2, -math.pi/3, math.pi/4])
    #tile_standalone.set_rotation([0, 0, 0]) 
    tile_standalone.set_remanence(1.2 / (4*math.pi*1e-7))
    tile_standalone.set_mag_angle([math.pi/3, math.pi/4])

    # Load reference points from COMSOL calculation
    path_COMSOL_eval = os.path.dirname(os.path.abspath(__file__)) + '/../util/evaluations/'

    with open(os.path.join(path_COMSOL_eval, 'Example_prism_normH_x_rot_rev.txt'), "r") as file:
                Tx = file.readlines()[8:]
    Tx_split = np.asarray([line.split() for line in Tx], dtype=np.float64)
    x = Tx_split[:,0]
    H_norm_x_COMSOL = Tx_split[:,1]
    off = np.ones(len(x))
    eval_points_x = np.c_[x, off*offset[1], off*offset[2]]

    # (updated_tiles_x, H_x) = MagTense.run_simulation(tile_val, eval_points_x)
    # Testing the standalone version
    (updated_tiles_x, H_dict_x) = MagTenseStandalone.run_simulation([tile_standalone], eval_points_x)
    H_x = np.array(list(H_dict_x.field.values()))

    H_norm_x_MagTense  = MagTense.get_norm_magnetic_flux(H_x)
    error = abs(H_norm_x_COMSOL - H_norm_x_MagTense)

    # Standard parameters in settings: max_error=0.00001, max_it=500
    # util_plot.create_plot(updated_tiles_x, eval_points_x, H_x)

    fig, ax = plt.subplots(1,3)
    fig.suptitle("MagTensePython vs. COMSOL")    
    ax[0].plot(x, H_norm_x_COMSOL, 'b', label='COMSOL')
    ax[0].plot(x, H_norm_x_MagTense, 'r', label='MagTense')
    # ax[0].plot(x, error, 'g', label= 'Error')
    ax[0].legend()
    ax[0].set_xlabel('x_axis')
    ax[0].set_ylabel('H_norm')


    with open(os.path.join(path_COMSOL_eval, 'Example_prism_normH_y_rot_rev.txt'), "r") as file:
                Ty = file.readlines()[8:]
    Ty_split = np.asarray([line.split() for line in Ty], dtype=np.float64)
    y = Ty_split[:,0]
    H_norm_y_COMSOL = Ty_split[:,1]
    off = np.ones(len(y))
    eval_points_y = np.c_[off*offset[0], y, off*offset[2]]

    # (updated_tiles_y, H_y) = MagTense.run_simulation(tile_val, eval_points_y)
    # Testing the standalone version
    (updated_tiles_y, H_dict_y) = MagTenseStandalone.run_simulation([tile_standalone], eval_points_y)
    H_y = np.array(list(H_dict_y.field.values()))

    H_norm_y_MagTense  = MagTense.get_norm_magnetic_flux(H_y)
    error = abs(H_norm_y_COMSOL - H_norm_y_MagTense)

    # Standard parameters in settings: max_error=0.00001, max_it=500
    # util_plot.create_plot(updated_tiles_y, eval_points_y, H_y)
  
    ax[1].plot(y, H_norm_y_COMSOL, 'b', label='COMSOL')
    ax[1].plot(y, H_norm_y_MagTense, 'r', label='MagTense')
    # ax[1].plot(y, error, 'g', label= 'Error')
    ax[1].legend()
    ax[1].set_xlabel('y_axis')
    ax[1].set_ylabel('H_norm')

    with open(os.path.join(path_COMSOL_eval, 'Example_prism_normH_z_rot_rev.txt'), "r") as file:
                Tz = file.readlines()[8:]
    Tz_split = np.asarray([line.split() for line in Tz], dtype=np.float64)
    z = Tz_split[:,0]
    H_norm_z_COMSOL = Tz_split[:,1]
    off = np.ones(len(z))
    eval_points_z = np.c_[off*offset[0], off*offset[1], z]

    # (updated_tiles_z, H_z) = MagTense.run_simulation(tile_val, eval_points_z)
    # Testing the standalone version
    (updated_tiles_z, H_dict_z) = MagTenseStandalone.run_simulation([tile_standalone], eval_points_z)
    H_z = np.array(list(H_dict_z.field.values()))

    H_norm_z_MagTense  = MagTense.get_norm_magnetic_flux(H_z)
    error = abs(H_norm_z_COMSOL - H_norm_z_MagTense)

    # Standard parameters in settings: max_error=0.00001, max_it=500
    # util_plot.create_plot(updated_tiles_z, eval_points_z, H_z)
  
    ax[2].plot(z, H_norm_z_COMSOL, 'b', label='COMSOL')
    ax[2].plot(z, H_norm_z_MagTense, 'r', label='MagTense')
    # ax[2].plot(z, error, 'g', label= 'Error')
    ax[2].legend()
    ax[2].set_xlabel('z_axis')
    ax[2].set_ylabel('H_norm')

    plt.show()

main()
