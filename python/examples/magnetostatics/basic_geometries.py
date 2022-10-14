#%%
import numpy as np

from magtense import magstatics
from magtense.utils.plot import create_plot


def basic_geometries():
    # Define all tiles
    all_tiles = magstatics.Tiles(6)
    all_tiles.M_rem= (1.2 / (4 * np.pi * 1e-7))
    all_tiles.tile_type = [2, 1, 3, 4, 5, 7]
    all_tiles.set_mag_angle()

    # 0: Prism
    all_tiles.size = ([0.1, 0.3, 0.2], 0)
    all_tiles.offset = ([0.1, 0.2, 0.1], 0)
    all_tiles.color = ([1, 0, 0], 0)

    # 1: Cylindrical Tiles
    all_tiles.center_pos = ([1, 0, 0.3], 1) 
    all_tiles.dev_center = ([0.15, np.pi/9, 0.3], 1)
    all_tiles.color = ([0, 0, 1], 1)

    # 2: Circpiece
    all_tiles.center_pos = ([0.85, np.pi/5, 1.2], 2) 
    all_tiles.dev_center = ([0.15, np.pi/7, 0.25], 2)
    all_tiles.color = ([1, 0.5, 0], 2)

    # 3: Inverted Circpiece
    all_tiles.center_pos = ([0.2, np.pi/6, 0.75], 3) 
    all_tiles.dev_center = ([0.05, np.pi/4, 0.4], 3)
    all_tiles.color = ([0.3, 0.8, 0.2], 3)

    # 4: Tetrahedron
    all_tiles.vertices = (np.array([[0.65, 0.9, 0.5],
                                    [0.8, 0.9, 0.7],
                                    [0.85, 0.55, 0.25],
                                    [0.95, 0.85, 0.15]]), 4)
    all_tiles.color = ([0, 0, 0], 4)

    # 5: Prolate Spheroid
    all_tiles.size = ([0.1, 0.3, 0.1], 5)
    all_tiles.offset = ([0.1, 0.6, 0.7], 5)
    all_tiles.rot = ([0, 0, 2], 5)
    all_tiles.color = ([1, 0, 1], 5)

    area = [0.9, 1, 0.8]
    n_points = [1, 100, 80]
    seg = area / np.asarray(n_points)
    pts = [[(i + 0.5) * seg[0], (j + 0.5) * seg[1], (k + 0.5) * seg[2]]
            for i in range(0,n_points[0])
            for j in range(0, n_points[1])
            for k in range(0, n_points[2])]
    pts = np.asarray(pts, dtype=np.float64, order='F')

    (updated_tiles, H_demag) = magstatics.run_simulation(all_tiles, pts)
    create_plot(updated_tiles, pts, H_demag)


if __name__ == '__main__':
    basic_geometries()

# %%
