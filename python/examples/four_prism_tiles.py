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
    places = [10, 10, 1]
    area = [1, 1, 0.01]
    # Defining occupied places in grid
    filled_positions = [[2, 4, 0], [4, 2, 0], [4, 6, 0], [6, 4, 0]]
    # Defining angle of magnetization in spherical coordinates (azimuth, polar angle) for each tile
    mag_angles = [[0, math.pi/2], [math.pi, math.pi/2], [math.pi, math.pi/2], [0, math.pi/2]]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = MagTense.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)
    print(tiles)

    # Standard parameters in settings: max_error=0.00001, max_it=500
    iterated_tiles = MagTense.iterate_magnetization(tiles)
    N = MagTense.get_N_tensor(iterated_tiles,points)
    H = MagTense.get_H_field(iterated_tiles,points,N)
    util_plot.create_plot(iterated_tiles, points, H, grid=grid)

    # All steps in one command - iterate_solution = True, return_field = True
    # N is not reused for calculating the H field
    # (updated_tiles, H) = MagTense.run_simulation(tiles, points, grid=grid, plot=False)
    # util_plot.create_plot(updated_tiles, points, H, grid=grid)

    # Testing the standalone version
    # (tiles, points, grid) = MagTenseStandalone.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)
    # (tiles, H) = MagTenseStandalone.run_simulation(tiles, points, grid=grid, plot=True)

main()
