import numpy as np

from magtense import magtense
from magtense.utils.eval import get_p2p, get_average_magnetic_flux


def halbach():
    # Defining grid
    places = [10, 10, 1]
    area = [1, 1, 0.01]

    # Defining filled places
    filled_positions = [[i, j, 0] for i in range(10) for j in range(10) \
        if (i==5 and j==5) or (i==4 and j==5) or (i==4 and j==4) or (i==5 and j==4)]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = magtense.setup(places, area, filled_positions=filled_positions)

    # Optimal angle - https://orbit.dtu.dk/ws/files/100219185/Analysis_of_the_magnetic_field_force_and_torque_for_two_dimensional_Halbach_cylinders.pdf
    for i in range(tiles.n):
        # Calculate angle of tile center to x-axis
        d = tiles.get_offset(i) - np.array([0.5, 0.5, 0.05])
        # r_comp = np.sqrt(d[0]**2, d[1]**2)
        angle_comp = np.arctan2(d[1], d[0])
        azimuth_opt = 2 * angle_comp
        tiles.set_mag_angle_i([np.pi / 2, azimuth_opt], i)

    # Standard parameters in settings: max_error=0.00001, max_it=500
    (updated_tiles, H) = magtense.run_simulation(tiles, points, grid=grid, plot=True)

    print("Average magnetic field: " + str(get_average_magnetic_flux(H)))
    print("Peak to peak: " + str(get_p2p(H)))


if __name__ == '__main__':
    halbach()