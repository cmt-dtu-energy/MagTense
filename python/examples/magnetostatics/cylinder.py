#%%
import math

from magtense import magstatics
from magtense.utils.plot import create_plot


def cylinder():
    spots = [3, 3, 3]
    area = [1,1,1]
    _, pts = magstatics.grid_config(spots, area, n_pts=[10,10,5])

    tile_val = magstatics.Tiles(
        n=1,
        center_pos=[0.3, math.pi/2, 0.5],
        dev_center=[0.3, math.pi/4, 0.1],
        offset=[0.8, -0.1, 0.3],
        rot=[0, 0, 0],
        # rot = [-math.pi/6, math.pi/5, math.pi/2],
        tile_type=1,
        M_rem=1.2/(4*math.pi*1e-7),
        color=[1, 0, 0],
        mag_angle=[math.pi/4, -math.pi/2]
    )

    it_tiles = magstatics.iterate_magnetization(tile_val)
    demag_tensor = magstatics.get_demag_tensor(it_tiles, pts)
    H_out = magstatics.get_H_field(it_tiles, pts, demag_tensor)
    
    create_plot(it_tiles, pts, H_out, spots, area)


if __name__ == '__main__':
    cylinder()

# %%
