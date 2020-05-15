import numpy as np

import MagTense

def main():
    # Define spherical tiles
    spheroids = MagTense.Tiles(3)

    # Defines first tile as sphere
    spheroids.set_tile_type_i(6, 0)
    spheroids.set_size_i([0.3, 0, 0], 0)
    spheroids.set_offset_i([1, 1, 1], 0)
    spheroids.set_remanence_i(1 / (4*np.pi*1e-7), 0)
    spheroids.set_easy_axis_i([1, 0, -1], 0)
    spheroids.set_color_i([1, 0, 0], 0)

    # Defines second tile as oblate spheroid
    spheroids.set_tile_type_i(7, 1)
    spheroids.set_size_i([0.1, 0.5, 0.5], 1)
    spheroids.set_offset_i([1, 2.5, 1], 1)
    spheroids.set_remanence_i(1 / (4*np.pi*1e-7), 1)
    spheroids.set_easy_axis_i([-1, 1, 0], 1)
    spheroids.set_rot_axis_i([-1, 1, -3], 'c', 1)
    spheroids.set_color_i([1, 0, 1], 1)

    # Defines thrid tile as prolate spheroid
    spheroids.set_tile_type_i(7, 2)
    spheroids.set_size_i([0.2, 0.5, 0.2], 2)
    spheroids.set_offset_i([2, 1, 1.5], 2)
    spheroids.set_remanence_i(1 / (4*np.pi*1e-7), 2)
    spheroids.set_easy_axis_i([1, 1, 1], 2)
    spheroids.set_rot_axis_i([1, 1, 1], 's', 2)
    spheroids.set_color_i([0, 0, 1], 2)

    # Defining the points, where the field is evaluated (output)
    area = [3, 3, 3]
    n_points = [7, 7, 7]
    seg = area/np.asarray(n_points)
    points = [[i*seg[0]+seg[0]/2, j*seg[1]+seg[1]/2, k*seg[2]+seg[2]/2] for i in range(0,n_points[0]) \
        for j in range(0, n_points[1]) for k in range(0, n_points[2])]
    points = np.asarray(points, dtype=np.float64, order='F')

    # Standard parameters in settings: max_error=0.00001, max_it=500
    (updated_tiles, H) = MagTense.run_simulation(spheroids, points, plot=True)

    print("Average magnetic field: " + str(MagTense.get_average_magnetic_flux(H)))
    print("Peak to peak: " + str(MagTense.get_p2p(H)))

if __name__ == '__main__':
    main()