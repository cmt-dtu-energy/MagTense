import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../lib_mag')
import numpy as np
import math

import MagTense
import MagTenseStandalone
import util_plot
import MagTenseSource

def main():
    # Defining grid and positions, angles of magnets
    places = [10, 10, 1]
    area = [1, 1, 0.01]
    filled_positions = [[2, 4, 0], [4, 2, 0], [4, 6, 0], [6, 4, 0]]
    mag_angles = [0, math.pi, math.pi, 0]

    # Standard parameters in settings: max_error=0.00001, max_it=500, iterate_solution=1, return_field=1
    # Optional parameters for setup: param_settings, n_magnets, area, eval_points, eval_mode 
    (tiles, points, grid) = MagTense.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)

    iterated_tiles = MagTense.iterate_magnetization(tiles)
    N = MagTense.get_N_tensor(iterated_tiles,points)
    H = MagTense.get_H_field(iterated_tiles,points,N)
    solution = MagTense.MagneticFieldIntensity(points, H)
    util_plot.create_plot(iterated_tiles, solution, grid=grid)

    # All steps in one command - iterate_solution = True, return_field = True
    # N is not reused for calculating the H field
    # (updated_tiles, solution) = MagTense.run_simulation(tiles, points, grid=grid, plot=False)
    # util_plot.create_plot(updated_tiles, solution, grid=grid)

    # Testing the standalone version
    # (tiles, points, grid) = MagTenseStandalone.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)
    # (tiles, H) = MagTenseStandalone.run_simulation(tiles, points, grid=grid, plot=True)

main()
