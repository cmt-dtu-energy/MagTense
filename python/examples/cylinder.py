import math

from magtense import magtense
from magtense.utils.plot import create_plot


def cylinder():
    B_rem = 1.2
    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = magtense.setup(places=[3, 3, 3], area=[1, 1, 1], eval_points=[10,10,5])

    # Define validation tile
    tile_val = magtense.Tiles(1)
    tile_val.set_center_pos([[0.3, math.pi/2, 0.5]])
    tile_val.set_dev_center([0.3, math.pi/4, 0.1])
    tile_val.set_offset_i([0.8, -0.1, 0.3],0)
    # tile_val.set_rotation_i([-math.pi/6, math.pi/5, math.pi/2],0)
    tile_val.set_rotation_i([0, 0, 0],0)
    tile_val.set_tile_type(1)
    tile_val.set_remanence(B_rem / (4*math.pi*1e-7))
    tile_val.set_mag_angle([[math.pi/4, -math.pi/2]])
    tile_val.set_color([[1, 0, 0]])

    # Standard parameters in settings: max_error=0.00001, max_it=500
    iterated_tiles = magtense.iterate_magnetization(tile_val)
    N = magtense.get_N_tensor(iterated_tiles,points)
    H = magtense.get_H_field(iterated_tiles,points,N)
    create_plot(iterated_tiles, points, H, grid=grid)


if __name__ == '__main__':
    cylinder()
