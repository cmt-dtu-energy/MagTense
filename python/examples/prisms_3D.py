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
    places = [10, 10, 5]
    area = [1, 1, 0.5]

    # Defining filled places
    filled_positions = [[4, 3, 4], [1, 2, 2], [3, 1, 0], [2, 7, 2], [8, 9, 3]]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = MagTense.setup(places, area, filled_positions=filled_positions, eval_points=[10, 10, 5])

    # Standard parameters in settings: max_error=0.00001, max_it=500
    (updated_tiles, H) = MagTense.run_simulation(tiles, points, grid=grid, plot=True)

    print("Average magnetic field: " + str(MagTense.get_average_magnetic_flux(H)))
    print("Peak to peak: " + str(MagTense.get_p2p(H)))

main()