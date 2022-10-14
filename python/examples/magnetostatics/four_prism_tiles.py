#%%
import math

from magtense import magstatics
from magtense.utils.plot import create_plot


def four_prisms():
    spots = [10, 10, 1]
    area = [1, 1, 0.01]
    filled_pos = [[2, 4, 0], [4, 2, 0], [4, 6, 0], [6, 4, 0]]

    # Defining angle of magnetization in spherical coordinates
    # (azimuth, polar angle) for each tile
    mag_angles = [[0, math.pi/2], [math.pi, math.pi/2],
                  [math.pi, math.pi/2], [0, math.pi/2]]

    tiles, pts = magstatics.grid_config(spots, area, filled_pos, mag_angles=mag_angles)
    print(tiles)

    it_tiles = magstatics.iterate_magnetization(tiles)
    demag_tensor = magstatics.get_demag_tensor(it_tiles, pts)
    H_out = magstatics.get_H_field(it_tiles, pts, demag_tensor)
    create_plot(it_tiles, pts, H_out, spots, area)


if __name__ == '__main__':
    four_prisms()

# %%
