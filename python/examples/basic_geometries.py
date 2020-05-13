import numpy as np

import MagTense

def main():
    # Define all tiles
    all_tiles = MagTense.Tiles(5)
    all_tiles.set_remanence(1.2 / (4*np.pi*1e-7))
    all_tiles.set_color([1, 0, 0])
    all_tiles.set_tile_type([2, 1, 3, 4, 5])
    all_tiles.set_mag_angle_rand()

    # 0: Prism
    all_tiles.set_size_i([0.1, 0.3, 0.2], 0) 
    all_tiles.set_offset_i([0.1, 0.5, 0.1], 0)

    # 1: Cylindrical Tiles
    all_tiles.set_center_pos_i([0.7, np.pi/8, 0.6], 1) 
    all_tiles.set_dev_center_i([0.15, np.pi/8, 0.4], 1)

    # 2: Circpiece
    all_tiles.set_center_pos_i([0.7, np.pi/3, 0.75], 2) 
    all_tiles.set_dev_center_i([0.15, np.pi/6, 0.3], 2)

    # 3: Inverted Circpiece
    all_tiles.set_center_pos_i([0.2, np.pi/6, 0.5], 3) 
    all_tiles.set_dev_center_i([0.05, np.pi/4, 0.25], 3)

    # 4: Tetrahedron
    all_tiles.set_vertices_i(np.array([[0.6,0.9,0.5],[0.7,0.8,0.7],[0.55,0.6,0.1],[0.9,0.9,0.15]]), 4)

    area = [1, 1, 1]
    n_points = [5, 5, 5]
    seg = area/np.asarray(n_points)
    points = [[i*seg[0]+seg[0]/2, j*seg[1]+seg[1]/2, k*seg[2]+seg[2]/2] for i in range(0,n_points[0]) \
        for j in range(0, n_points[1]) for k in range(0, n_points[2])]
    points = np.asarray(points, dtype=np.float64, order='F')

    # Simulation
    (updated_tiles, H) = MagTense.run_simulation(all_tiles, points, plot=True)

if __name__ == '__main__'
    main()
