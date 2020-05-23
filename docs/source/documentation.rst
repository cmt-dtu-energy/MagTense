Documentation
========================================

========================================
Python
========================================

----------------------------------------
Objects
----------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`MagTiles <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/source/MagTense.py#L9>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    def __init__(self, n):
        # Initialization of arrays for specific tile parameters
        # Input to Fortran derived type MagTile     
        self.center_pos = np.zeros(shape=(n,3), dtype=np.float64, order='F') # r0, theta0, z0
        self.dev_center = np.zeros(shape=(n,3), dtype=np.float64, order='F') # dr, dtheta, dz
        self.size = np.zeros(shape=(n,3), dtype=np.float64, order='F') # a, b, c
        self.vertices = np.zeros(shape=(n,3,4), dtype=np.float64, order='F') # v1, v2, v3, v4 as column vectors
        self.M = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Mx, My, Mz
        self.u_ea = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Easy axis
        self.u_oa1 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.u_oa2 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.mu_r_ea = np.ones(shape=(n), dtype=np.float64, order='F')
        self.mu_r_oa = np.ones(shape=(n), dtype=np.float64, order='F')
        self.M_rem = np.zeros(shape=(n), dtype=np.float64, order='F')
        # Implemented tile types: 
        # 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv,
        # 5 = tetrahedron, 6 = sphere, 7 = spheroid, 10 = ellipsoid
        self.tile_type = np.ones(n, dtype=np.int32, order='F')
        self.offset = np.zeros(shape=(n,3), dtype=np.float64, order='F') # offset of global coordinates
        self.rot = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.color = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.magnetic_type = np.ones(n, dtype=np.int32, order='F') # 1 = hard magnet, 2 = soft magnet
        self.stfcn_index = np.ones(shape=(n), dtype=np.int32, order='F') # default index into the state function
        self.incl_it = np.ones(shape=(n), dtype=np.int32, order='F') # if equal to zero the tile is not included in the iteration
        self.use_sym = np.zeros(shape=(n), dtype=np.int32, order='F') # whether to exploit symmetry
        self.sym_op = np.ones(shape=(n,3), dtype=np.float64, order='F') # 1 for symmetry and -1 for anti-symmetry respectively to the planes
        self.M_rel = np.zeros(shape=(n), dtype=np.float64, order='F')
 
        # Internal parameters for python to prepare configuration
        self.grid_pos = np.zeros(shape=(n,3), dtype=np.float64, order='F') # positions in the grid
        self.n = n

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Evaluation points
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    points = np.zeros(shape=(number_of_points,3), dtype=np.float64, order='F')

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Demagnetization tensor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    N = np.zeros(shape=(number_of_tiles,number_of_points,3,3), dtype=np.float64, order='F')

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
H-field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    H = np.zeros(shape=(number_of_points,3), dtype=np.float64, order='F')


----------------------------------------
Functions
----------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`iterate_magnetization (tiles, **options) <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/source/MagTense.py#L441>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Returns updated tiles. Iterates through the given tiles to determine their influence 
on each other.

::

    # Options:
    max_error = 0.00001 # Iteration stops if magnetization change of tiles is below this value
    max_it = 500 # Maximum number of performed iterations
    T = 300. # Temperature for the state function if required

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`get_N_tensor (tiles, points) <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/source/MagTense.py#L459>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Returns the demagnetization tensor N of the given tiles and 
the specified evaluation points.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`get_H_field (tiles, points, N=None) <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/source/MagTense.py#L468>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Returns the magnetic field H at the specified evaluation 
points of the given tiles.
Optionally, a precalculated demagnetization tensor N can be 
handed over in order to prevent unnecessary and expensive 
recalculation of N if geometry of the setup does not change.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`run_simulation (tiles, points, **options) <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/source/MagTense.py#L406>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Does the previous mentioned steps all together.
Demagnetization tensor will not be reused.

::

    # Options:
    plot = False # Boolean if results shall be plotted
    grid = None # Optional grid can be displayed
    max_error = 0.00001 # Iteration stops if magnetization change of tiles is below this value
    max_it = 500 # Maximum number of performed iterations
    T = 300. # Temperature for the state function if required
    iterate_solution = True # Boolean if the magnetization of the tiles shall be iterated
    return_field=True # Boolean if magnetic field shall be calculated

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`create_plot (iterated_tiles, points, H, grid=None) <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/source/util_plot.py#L422>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Creates a matplotlib with the iterated tiles and the calculated 
H-field at the evaluation points as quiver plot.
Additionally, an optional grid can be displayed.