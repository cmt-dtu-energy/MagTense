#%%
import numpy as np

from magtense import magstatics
from magtense.utils.eval import get_p2p, get_average_magnetic_flux
from magtense.utils.plot import create_plot


def prism_grid():
    places = [10, 10, 5]
    area = [1, 1, 0.5]
    filled_pos = [[4, 3, 4], [1, 2, 2], [3, 1, 0], [2, 7, 2], [8, 9, 3]]

    tiles, points = magstatics.grid_config(places, area, filled_pos, [10, 10, 5])

    # Standard parameters in settings: max_error=0.00001, max_it=500
    (_, H_out) = magstatics.run_simulation(tiles, points)

    print("Average magnetic field: " + str(get_average_magnetic_flux(H_out)))
    print("Peak to peak: " + str(get_p2p(H_out)))


def prism_multiple(n_mag=1, soft=None, res=16, x_max=1, y_max=1, z_max=1):
    mu0 = 4 * np.pi * 1e-7
    a = 0.1
    b = 0.3
    c = 0.2
    
    if soft is None: soft = [0 for _ in range(n_mag)]

    prism = magstatics.Tiles(n_mag, tile_type=2)

    for i in range(n_mag):
        # Relative permeability
        if soft == 1:
            prism.mu_r_ea = (4000, i)
            prism.mu_r_oa = (4000, i)
            prism.remanence = (0, i)
            prism.mag_type = (2, i)

            # Default: Easy axis in direction of x-axis
            # [polar angle, azimuthal angle]
            prism.set_mag_angle_i([np.pi/2, 0], i)

        else:
            prism.mu_r_ea = (1.06, i)
            prism.mu_r_oa = (1.17, i)

            # Set remanence in [A/m]
            prism.M_rem = (1.2 / mu0, i)

            # Default: Easy axis in direction of x-axis
            # [polar angle, azimuthal angle]
            prism.set_mag_angle_i([np.random.rand() * np.pi, np.random.rand() * 2 * np.pi], i)

        prism.size = ([a, b, c], i)
        prism.offset = ([0.1 + i*0.5, 0.2 + i*0.5, 0.1 + i*0.5], i)
        prism.set_mag_angle()
        prism.color = ([i / n_mag, 0, 0], i)

    x_eval = np.linspace(0, x_max, res)
    y_eval = np.linspace(0, y_max, res)
    z_eval = np.linspace(0, z_max, res)
    xv, yv, zv = np.meshgrid(x_eval, y_eval, z_eval)
    pts_eval = np.hstack([xv.reshape(-1,1), yv.reshape(-1,1), zv.reshape(-1,1)])

    (updated_tiles, H_out) = magstatics.run_simulation(prism, pts_eval)
    create_plot(updated_tiles, pts_eval, H_out)


if __name__ == '__main__':
    # prism_grid()
    prism_multiple(n_mag=2, res=8, soft=[0,1])
# %%
