import numpy as np

from magtense import magtense
from magtense.utils.plot import create_plot


def basic_geometries():
    # Define all tiles
    all_tiles = magtense.Tiles(6)
    all_tiles.set_remanence(1.2 / (4 * np.pi * 1e-7))
    all_tiles.set_tile_type([2, 1, 3, 4, 5, 7])
    all_tiles.set_mag_angle_rand()

    # 0: Prism
    all_tiles.set_size_i([0.1, 0.3, 0.2], 0) 
    all_tiles.set_offset_i([0.1, 0.2, 0.1], 0)
    all_tiles.set_color_i([1, 0, 0], 0)

    # 1: Cylindrical Tiles
    all_tiles.set_center_pos_i([1, 0, 0.3], 1) 
    all_tiles.set_dev_center_i([0.15, np.pi/9, 0.3], 1)
    all_tiles.set_color_i([0, 0, 1], 1)

    # 2: Circpiece
    all_tiles.set_center_pos_i([0.85, np.pi/5, 1.2], 2) 
    all_tiles.set_dev_center_i([0.15, np.pi/7, 0.25], 2)
    all_tiles.set_color_i([1, 0.5, 0], 2)

    # 3: Inverted Circpiece
    all_tiles.set_center_pos_i([0.2, np.pi/6, 0.75], 3) 
    all_tiles.set_dev_center_i([0.05, np.pi/4, 0.4], 3)
    all_tiles.set_color_i([0.3, 0.8, 0.2], 3)

    # 4: Tetrahedron
    all_tiles.set_vertices_i(np.array([[0.65, 0.9, 0.5],
                                       [0.8, 0.9, 0.7],
                                       [0.85, 0.55, 0.25],
                                       [0.95, 0.85, 0.15]]), 4)
    all_tiles.set_color_i([0, 0, 0], 4)

    # 5: Prolate Spheroid
    all_tiles.set_size_i([0.1, 0.3, 0.1], 5)
    all_tiles.set_offset_i([0.1, 0.6, 0.7], 5)
    all_tiles.set_rot_axis_i([0, 0, 2], 5)
    all_tiles.set_color_i([1, 0, 1], 5)

    area = [0.9, 1, 0.8]
    n_points = [1, 100, 80]
    seg = area / np.asarray(n_points)
    points = [[i * seg[0] + seg[0] / 2,
               j * seg[1] + seg[1] / 2,
               k * seg[2] + seg[2] / 2]
               for i in range(0,n_points[0])
               for j in range(0, n_points[1])
               for k in range(0, n_points[2])]
    points = np.asarray(points, dtype=np.float64, order='F')

    # Simulation
    (updated_tiles, H) = magtense.run_simulation(all_tiles, points, plot=False)

    # Plotting
    create_plot(updated_tiles, H=H, eval_points=points)


if __name__ == '__main__':
    basic_geometries()
