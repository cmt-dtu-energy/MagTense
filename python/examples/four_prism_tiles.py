import math

from magtense import magtense
from magtense.utils.plot import create_plot


def four_prisms():
    # Defining grid
    places = [10, 10, 1]
    area = [1, 1, 0.01]
    # Defining occupied places in grid
    filled_positions = [[2, 4, 0], [4, 2, 0], [4, 6, 0], [6, 4, 0]]
    # Defining angle of magnetization in spherical coordinates (azimuth, polar angle) for each tile
    mag_angles = [[0, math.pi/2], [math.pi, math.pi/2], [math.pi, math.pi/2], [0, math.pi/2]]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = magtense.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)
    print(tiles)

    # Standard parameters in settings: max_error=0.00001, max_it=500
    iterated_tiles = magtense.iterate_magnetization(tiles)
    N = magtense.get_N_tensor(iterated_tiles,points)
    H = magtense.get_H_field(iterated_tiles,points,N)
    create_plot(iterated_tiles, points, H, grid=grid)

    # All steps in one command - iterate_solution = True, return_field = True
    # N is not reused for calculating the H field
    # (updated_tiles, H) = magtense.run_simulation(tiles, points, grid=grid, plot=False)
    # create_plot(updated_tiles, points, H, grid=grid)


if __name__ == '__main__':
    four_prisms()
