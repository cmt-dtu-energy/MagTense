import numpy as np

import magtense
from util_plot import create_plot

def prism_grid():
    # Defining grid
    places = [10, 10, 5]
    area = [1, 1, 0.5]

    # Defining filled places
    filled_positions = [[4, 3, 4], [1, 2, 2], [3, 1, 0], [2, 7, 2], [8, 9, 3]]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = magtense.setup(places, area, filled_positions=filled_positions, eval_points=[10, 10, 5])

    # Standard parameters in settings: max_error=0.00001, max_it=500
    (updated_tiles, H) = magtense.run_simulation(tiles, points, grid=grid, plot=True)

    print("Average magnetic field: " + str(magtense.get_average_magnetic_flux(H)))
    print("Peak to peak: " + str(magtense.get_p2p(H)))


def prism_multiple(n_mag=1, soft=None, res=16, x_max=1, y_max=1, z_max=1):
    mu0 = 4 * np.pi * 1e-7
    a = 0.1
    b = 0.3
    c = 0.2
    
    if soft is None: soft = [0 for _ in range(n_mag)]

    prism = magtense.Tiles(n_mag)
    prism.set_tile_type(2)

    for i in range(n_mag):
        # Relative permeability
        if soft == 1:
            prism.set_mu_r_ea_i(4000, i)
            prism.set_mu_r_oa_i(4000, i)
            prism.set_remanence_i(0, i)
            prism.set_mag_type_i(2, i)

            # Default: Easy axis in direction of x-axis
            # [polar angle, azimuthal angle]
            prism.set_mag_angle_i([np.pi/2, 0], i)

        else:
            prism.set_mu_r_ea_i(1.06, i)
            prism.set_mu_r_oa_i(1.17, i)

            # Set remanence in [A/m]
            prism.set_remanence_i(1.2 / mu0, i)

            # Default: Easy axis in direction of x-axis
            # [polar angle, azimuthal angle]
            prism.set_mag_angle_i([np.random.rand() * np.pi, np.random.rand() * 2 * np.pi], i)

        prism.set_size_i([a, b, c], i) 
        prism.set_offset_i([0.1 + i*0.5, 0.2 + i*0.5, 0.1 + i*0.5], i)
        prism.set_mag_angle_rand()
        prism.set_color_i([i / n_mag, 0, 0], i)

    x_eval = np.linspace(0, x_max, res)
    y_eval = np.linspace(0, y_max, res)
    z_eval = np.linspace(0, z_max, res)
    xv, yv, zv = np.meshgrid(x_eval, y_eval, z_eval)
    pts_eval = np.hstack([xv.reshape(-1,1), yv.reshape(-1,1), zv.reshape(-1,1)])

    # Simulation
    (updated_tiles, H) = magtense.run_simulation(prism, pts_eval, plot=False)

    # Plotting
    create_plot(updated_tiles, H=H, eval_points=pts_eval)


if __name__ == '__main__':
    # prism_grid()
    prism_multiple(n_mag=2, res=8, soft=[0,1])