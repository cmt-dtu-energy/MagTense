# Python Interface of MagTense

The fortran code can be directly called from Python.
The tool f2py of the numpy package was used to wrap the interface file *lib_mag/FortranToPythonIO.f90*.

Accessible functions:

- **iterate_magnetization**(tiles, **options):  
    Iterates through the given tiles to determine their influence on each other.
    Updated tiles are returned.

- **get_N_tensor**(tiles, points):  
    Returns the demagnetization tensor N of the given tiles and the specified evaluation points.

- **get_H_field**(tiles, points, N=None):  
    Returns the magentic field H at the specified evaluation points of the given tiles.  
    Optionally, a precalculated demagnization tensor N can be handed over in order to prevent unnecessary and expensive recalculation of N iuf geometry of the setup does not change.

- **create_plot**(iterated_tiles, H, grid=None):  
    Creates a matlibplot with the iterated tiles and the calculated magnetic field H as quiver plot.  
    Additionally, an optional grid can be displayed.

- **run_simulation**(tiles, points, grid=None, plot=False, **options):  
    Does all the previous mentioned steps all together.  
    Demagnetization tensor will not be reused.

The required packages to run the code can be found in *documentation/requirements.txt*.  
In order to run the code the Fortran code has to be wrapped with f2py on your machine. The usage of the corresponding Makefile can be found in subfolder *lib_mag/*.
