# Python Interface of MagTense

The fortran code can be directly called from Python.
The tool f2py of the numpy package is used to wrap the interface file *lib_mag/FortranToPythonIO.f90*.

## Accessible functions

- **iterate_magnetization**(tiles, **options):  
    Iterates through the given tiles to determine their influence on each other.  
    Updated tiles are returned.

    ```Python
    # Options:
    max_error = 0.00001 # Iteration stops if magnetization change of tiles is below this value
    max_it = 500 # Maximum number of performed iterations
    T = 300. # Temperature for the state function if required
    ```

- **get_N_tensor**(tiles, points):  
    Returns the demagnetization tensor N of the given tiles and the specified evaluation points.

- **get_H_field**(tiles, points, N=None):  
    Returns the magentic field H at the specified evaluation points of the given tiles.  
    Optionally, a precalculated demagnization tensor N can be handed over in order to prevent unnecessary and expensive recalculation of N if geometry of the setup does not change.

- **create_plot**(iterated_tiles, points, H, grid=None):  
    Creates a matlibplot with the iterated tiles and the calculated magnetic field H at the evaluation points as quiver plot.  
    Additionally, an optional grid can be displayed.

- **run_simulation**(tiles, points, **options):  
    Does all the previous mentioned steps all together.  
    Demagnetization tensor will not be reused.

    ```Python
    # Options:
    plot = False # Boolean if results shall be plotted
    grid = None # Optional grid can be displayed
    max_error = 0.00001 # Iteration stops if magnetization change of tiles is below this value
    max_it = 500 # Maximum number of performed iterations
    T = 300. # Temperature for the state function if required
    iterate_solution = True # Boolean if the magnetization of the tiles shall be iterated
    return_field=True # Boolean if magnetic field shall be calculated
    ```

## Arguments

- Magnetic Tiles [Defined in **source/MagTense.py** as **Tiles()**]

    ```Python
    def __init__(self, n):
        # Initialization of arrays for specific tile parameters
        # Input to Fortran derived type MagTile
        self.center_pos = np.zeros(shape=(n,3), dtype=np.float64, order='F') # r0, theta0, z0
        self.dev_center = np.zeros(shape=(n,3), dtype=np.float64, order='F') # dr, dtheta, dz
        self.size = np.zeros(shape=(n,3), dtype=np.float64, order='F') # a, b, c
        self.M = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Mx, My, Mz
        self.u_ea = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Easy axis
        self.u_oa1 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.u_oa2 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.mu_r_ea = np.ones(shape=(n), dtype=np.float64, order='F')
        self.mu_r_oa = np.ones(shape=(n), dtype=np.float64, order='F')
        self.M_rem = np.zeros(shape=(n), dtype=np.float64, order='F')
        self.tile_type = np.ones(n, dtype=np.int32, order='F') # 1 = cylinder, 2 = prism, 3 = ellipsoid
        self.offset = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Offset of global coordinates
        self.rot = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.color = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self.magnetic_type = np.ones(n, dtype=np.int32, order='F') # 1 = hard magnet, 2 = soft magnet
        self.stfcn_index = np.ones(shape=(n), dtype=np.int32, order='F') # Default index into the state function
        self.incl_it = np.ones(shape=(n), dtype=np.int32, order='F') # If equal to zero the tile is not included in the iteration
        self.use_sym = np.ones(shape=(n), dtype=np.int32, order='F') # Whether to exploit symmetry
        self.sym_op = np.ones(shape=(n,3), dtype=np.float64, order='F') # 1 for symmetry and -1 for anti-symmetry respectively to the planes
        self.M_rel = np.zeros(shape=(n), dtype=np.float64, order='F')
        self.grid_pos = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Positions in the grid
        self.n = n
    ```

- Evaluation Points

    ```Python
        points = np.zeros(shape=(number_of_points,3), dtype=np.float64, order='F')
    ```

- Demagnetization Tensor

    ```Python
        N = np.zeros(shape=(number_of_tiles,number_of_points,3,3), dtype=np.float64, order='F')
    ```

- Magnetic Field

    ```Python
        H = np.zeros(shape=(number_of_points,3), dtype=np.float64, order='F')
    ```

- Grid [Defined in **source/MagTense.py** as **Grid()**]

## Additional remarks

The required packages to run the code can be found in *documentation/requirements.txt*.  
In order to run the code the Fortran code has to be wrapped with f2py on your machine. The usage of the corresponding Makefile can be found in subfolder *lib_mag/*.
