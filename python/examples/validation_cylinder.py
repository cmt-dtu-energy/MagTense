import os
import sys
import numpy as np
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../source')
import MagTense
import MagTenseStandalone
import util_plot

def main():
    # Defining grid
    places = [3, 3, 3]
    area = [1, 1, 1]
    B_rem = 1.2

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = MagTense.setup(places, area, eval_points=[10,10,5])

    # Define validation tile
    tile_val = MagTense.Tiles(1)
    tile_val.set_center_pos([[0.3, math.pi/2, 0.25]])
    tile_val.set_dev_center([0.2, math.pi/4, 0.4])
    tile_val.set_offset_i([0.5, 0.5, 0.5],0)
    tile_val.set_rotation_i([0, 0, 0],0)
    tile_val.set_tile_type(1)
    tile_val.set_remanence(B_rem / (4*math.pi*1e-7))
    tile_val.set_mag_angle([[math.pi/8, math.pi/2]])
    tile_val.set_color([[1, 0, 0]])

    # Standard parameters in settings: max_error=0.00001, max_it=500
    iterated_tiles = MagTense.iterate_magnetization(tile_val)
    N = MagTense.get_N_tensor(iterated_tiles,points)
    H = MagTense.get_H_field(iterated_tiles,points,N)
    util_plot.create_plot(None, points, H, grid=grid)

main()
