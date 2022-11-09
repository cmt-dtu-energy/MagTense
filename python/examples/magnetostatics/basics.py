#%%
import numpy as np

from magtense import magstatics
from magtense.utils import create_plot

mu0 = 4 * np.pi * 1e-7

def basic_geometries():
    tiles = magstatics.Tiles(
        n=6,
        M_rem=1.2 / mu0,
        tile_type=[2, 1, 3, 4, 5, 7],
        color=[[1, 0, 0], [0, 0, 1], [1, 0.5, 0],
               [0.3, 0.8, 0.2], [0, 0, 0], [1, 0, 1]]
    )

    # 0: Prism
    tiles.size = ([0.1, 0.3, 0.2], 0)
    tiles.offset = ([0.1, 0.2, 0.1], 0)

    # 1: Cylindrical Tiles
    tiles.center_pos = ([1, 0, 0.3], 1) 
    tiles.dev_center = ([0.15, np.pi/9, 0.3], 1)

    # 2: Circpiece
    tiles.center_pos = ([0.85, np.pi/5, 1.2], 2) 
    tiles.dev_center = ([0.15, np.pi/7, 0.25], 2)

    # 3: Inverted Circpiece
    tiles.center_pos = ([0.2, np.pi/6, 0.75], 3) 
    tiles.dev_center = ([0.05, np.pi/4, 0.4], 3)

    # 4: Tetrahedron
    tiles.vertices = (np.array([[0.65, 0.9, 0.5],
                                [0.8, 0.9, 0.7],
                                [0.85, 0.55, 0.25],
                                [0.95, 0.85, 0.15]]), 4)

    # 5: Prolate Spheroid
    tiles.size = ([0.1, 0.3, 0.1], 5)
    tiles.offset = ([0.1, 0.6, 0.7], 5)
    tiles.rot = ([0, 0, 2], 5)

    n_points = [1, 100, 80]
    seg = [0.9, 1, 0.8] / np.asarray(n_points)
    pts = [[(i + 0.5) * seg[0], (j + 0.5) * seg[1], (k + 0.5) * seg[2]]
            for i in range(0, n_points[0])
            for j in range(0, n_points[1])
            for k in range(0, n_points[2])]
    pts = np.asarray(pts, dtype=np.float64, order='F')

    updated_tiles, H_demag = magstatics.run_simulation(tiles, pts)
    create_plot(updated_tiles, pts, H_demag)


def cylinder():
    spots = [3, 3, 3]
    area = [1, 1, 1]
    _, pts = magstatics.grid_config(spots, area, n_pts=[10,10,5])

    tile_val = magstatics.Tiles(
        n=1,
        center_pos=[0.3, np.pi/2, 0.5],
        dev_center=[0.3, np.pi/4, 0.1],
        offset=[0.8, -0.1, 0.3],
        rot=[0, 0, 0],
        # rot=[-np.pi/6, np.pi/5, np.pi/2],
        tile_type=1,
        M_rem=1.2 / mu0,
        color=[1, 0, 0],
        mag_angle=[np.pi/4, -np.pi/2]
    )

    it_tiles = magstatics.iterate_magnetization(tile_val)
    demag_tensor = magstatics.get_demag_tensor(it_tiles, pts)
    H_out = magstatics.get_H_field(it_tiles, pts, demag_tensor)    
    create_plot(it_tiles, pts, H_out, spots, area)


def spheroid():
    spheroids = magstatics.Tiles(
        n=3,
        tile_type=[6,7,7],
        size=[[0.3, 0, 0], [0.1, 0.5, 0.5], [0.2, 0.5, 0.2]],
        offset=[[1, 1, 1], [1, 2.5, 1], [2, 1, 1.5]],
        M_rem=1 / mu0,
        easy_axis=[[1, 0, -1], [-1, 1, 0], [1, 1, 1]],
        rot=[[0, 0, 0], [-1, 1, -3], [1, 1, 1]],
        color=[[1, 0, 0], [1, 0, 1], [0, 0, 1]]
    )

    n_points = [7, 7, 7]
    seg = [3, 3, 3] / np.asarray(n_points)
    points = [[i * seg[0] + seg[0]/2,
               j * seg[1] + seg[1]/2,
               k * seg[2] + seg[2]/2]
               for i in range(0, n_points[0])
               for j in range(0, n_points[1])
               for k in range(0, n_points[2])]
    points = np.asarray(points, dtype=np.float64, order='F')

    updated_tiles, H_out = magstatics.run_simulation(spheroids, points)
    create_plot(updated_tiles, points, H_out)

    norm = [np.linalg.norm(H_point) * mu0 for H_point in H_out]

    print(f'Average magnetic field: {sum(norm) / len(norm):.3f}')
    print(f'Peak to peak: {max(norm) - min(norm):.3f}')


def four_prisms():
    spots = [10, 10, 1]
    area = [1, 1, 0.01]

    tiles, pts = magstatics.grid_config(
        spots,
        area,
        filled_pos=[[2, 4, 0], [4, 2, 0],
                    [4, 6, 0], [6, 4, 0]],
        mag_angles=[[0, np.pi/2], [np.pi, np.pi/2],
                    [np.pi, np.pi/2], [0, np.pi/2]]
    )
    print(tiles)

    it_tiles = magstatics.iterate_magnetization(tiles)
    demag_tensor = magstatics.get_demag_tensor(it_tiles, pts)
    H_out = magstatics.get_H_field(it_tiles, pts, demag_tensor)
    create_plot(it_tiles, pts, H_out, spots, area)


def prism_grid():
    tiles, points = magstatics.grid_config(
        spots=[10, 10, 5],
        area=[1, 1, 0.5],
        filled_pos=[[4, 3, 4], [1, 2, 2], [3, 1, 0],
                    [2, 7, 2], [8, 9, 3]],
        n_pts=[10, 10, 5]
    )
    _, H_out = magstatics.run_simulation(tiles, points)
    norm = [np.linalg.norm(H_point) * mu0 for H_point in H_out]

    print(f'Average magnetic field: {sum(norm) / len(norm):.3f}')
    print(f'Peak to peak: {max(norm) - min(norm):.3f}')


def prism_multiple(n_mag=1, soft=None, res=16, max=[1,1,1]):
    if soft is None: soft = [0 for _ in range(n_mag)]
    prism = magstatics.Tiles(n_mag, tile_type=2, size=[0.1, 0.3, 0.2])

    for i in range(n_mag):
        if soft == 1:
            prism.mu_r_ea = (4000, i)
            prism.mu_r_oa = (4000, i)
            prism.M_rem = (0, i)
            prism.magnet_type = (2, i)
            prism.set_easy_axis([np.pi/2, 0], i)

        else:
            prism.mu_r_ea = (1.06, i)
            prism.mu_r_oa = (1.17, i)
            prism.M_rem = (1.2 / mu0, i)
            prism.set_easy_axis(idx=i)

        prism.offset = ([0.1 + i/2, 0.2 + i/2, 0.1 + i/2], i)
        prism.color = ([i / n_mag, 0, 0], i)

    l_n = complex(0, res)
    coords = np.mgrid[0:max[0]:l_n, 0:max[1]:l_n, 0:max[2]:l_n]
    # Transpose from 'ij' to 'xy' notation
    pts_eval = coords.transpose((0,2,1,3)).reshape(3,-1).T

    updated_tiles, H_out = magstatics.run_simulation(prism, pts_eval)
    create_plot(updated_tiles, pts_eval, H_out)


def halbach_prism():
    '''
    Magnetic tiles shaped as rectangular prisms in Halbach configuration.
    Optimal angle:
    https://orbit.dtu.dk/ws/files/100219185/Analysis_of_the_magnetic_field_force_and_torque_for_two_dimensional_Halbach_cylinders.pdf
    '''
    filled_pos = [[i, j, 0] for i in range(10) for j in range(10) \
        if (i==5 and j==5) or (i==4 and j==5) or (i==4 and j==4) or (i==5 and j==4)]
    tiles, points = magstatics.grid_config([10, 10, 1], [1, 1, 0.01], filled_pos)

    for i in range(tiles.n):
        d = tiles.offset[i] - np.array([0.5, 0.5, 0.05])
        # r_comp = np.sqrt(d[0]**2, d[1]**2)
        angle_comp = np.arctan2(d[1], d[0])
        azimuth_opt = 2 * angle_comp
        tiles.set_easy_axis([np.pi / 2, azimuth_opt], i)

    _, H_out = magstatics.run_simulation(tiles, points)
    norm = [np.linalg.norm(H_point) * mu0 for H_point in H_out]
    print(f'Average magnetic field: {sum(norm) / len(norm):.3f}')
    print(f'Peak to peak: {max(norm) - min(norm):.3f}')
# %%
