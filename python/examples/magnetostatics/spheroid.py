#%%
import numpy as np

from magtense import magstatics
from magtense.utils.eval import get_p2p, get_average_magnetic_flux
from magtense.utils.plot import create_plot


def spheroid():
    # Define spherical tiles
    spheroids = magstatics.Tiles(3)

    # Defines first tile as sphere
    spheroids.tile_type = (6, 0)
    spheroids.size = ([0.3, 0, 0], 0)
    spheroids.offset = ([1, 1, 1], 0)
    spheroids.M_rem = (1 / (4*np.pi*1e-7), 0)
    spheroids.u_ea = ([1, 0, -1], 0)
    spheroids.color = ([1, 0, 0], 0)

    # Defines second tile as oblate spheroid
    spheroids.tile_type = (7, 1)
    spheroids.size = ([0.1, 0.5, 0.5], 1)
    spheroids.offset = ([1, 2.5, 1], 1)
    spheroids.M_rem = (1 / (4*np.pi*1e-7), 1)
    spheroids.u_ea = ([-1, 1, 0], 1)
    spheroids.rot = ([-1, 1, -3], 1)
    spheroids.color = ([1, 0, 1], 1)

    # Defines thrid tile as prolate spheroid
    spheroids.tile_type = (7, 2)
    spheroids.size = ([0.2, 0.5, 0.2], 2)
    spheroids.offset = ([2, 1, 1.5], 2)
    spheroids.M_rem = (1 / (4*np.pi*1e-7), 2)
    spheroids.u_ea = ([1, 1, 1], 2)
    spheroids.rot = ([1, 1, 1], 2)
    spheroids.color = ([0, 0, 1], 2)

    # Defining the points, where the field is evaluated (output)
    area = [3, 3, 3]
    n_points = [7, 7, 7]
    seg = area / np.asarray(n_points)
    points = [[i * seg[0] + seg[0]/2,
               j * seg[1] + seg[1]/2,
               k * seg[2] + seg[2]/2]
               for i in range(0,n_points[0])
               for j in range(0, n_points[1])
               for k in range(0, n_points[2])]
    points = np.asarray(points, dtype=np.float64, order='F')

    (updated_tiles, H_out) = magstatics.run_simulation(spheroids, points)
    create_plot(updated_tiles, points, H_out)

    print("Average magnetic field: " + str(get_average_magnetic_flux(H_out)))
    print("Peak to peak: " + str(get_p2p(H_out)))


if __name__ == '__main__':
    spheroid()
# %%
