<div align="center">
  <img src="https://cmt-dtu-energy.github.io/MagTense/_static/MagTense_logo.PNG" height=250><br>
</div>

# MagTense version 2021

MagTense consists of both a magnetostatic and a micromagnetism calculation framework.

The magnetostatic framework can calculate the magnetic field from objects shaped as cylinders, pieces of cylinders, prisms, circular pieces and tetrahedrons. This is done using a fully analytical calculation of the demagnetization tensor. The framework is fully implemented in Fortran and has both a Matlab MEX interface and a Python interface.

The micromagnetism framework solves the Landau-Lifshitz equation. The framework is fully implemented in Fortran and has a Matlab MEX interface and a Python interface, as well as an older Matlab implementation. The micromagnetism framework utilizes the magnetostatic framework for calculating the demagnetization field.

The webpage of the code is available at https://www.magtense.org.

The TechManual on the code is available at https://cmt-dtu-energy.github.io/MagTense/.


## Usage with Matlab

MagTense is directly useable in Matlab on Windows by downloading the already compiled MEX-files in [matlab/MEX](https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/MEX_files). The files are directly useable with no compilation required, although Matlab 2018a or greater is required. Examples of how to calculate magnetostatic and micromagnetic problems using the Matlab interface can be found in [matlab/examples](https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples).


### Compilation with a Visual Studio project file

If you want to compile MagTense with a Visual Studio project file for Windows, [MagTense.sln](https://github.com/cmt-dtu-energy/MagTense/blob/master/MagTense.sln), is available, as well as a Matlab function to build the MEX-files, [buildMagTenseMEX.m](https://github.com/cmt-dtu-energy/MagTense/blob/master/buildMagTenseMEX.m). MagTense utilizes Intel MKL for the micromagnetic simlations and can also utilize CUDA and CVODE.


## Usage with Python interface

Instructions on how to use the Python interface are listed in [python/](https://github.com/cmt-dtu-energy/MagTense/tree/master/python). Binary installers for the Python interface are available at the [Python
Package Index (PyPI)](https://pypi.org/project/magtense) and can be installed via (requires **Python 3.9** or **3.10**):

```sh
pip install magtense
```

### [Linux] Further requirements
- Install MKL and runtime for Intel® Fortran Compiler
    ```sh
    conda install mkl intel-fortran-rt
    ```
- Set `LD_PRELOAD` to selected mkl libraries and `LD_LIBRARY_PATH` to `lib` of active conda environment. A guide on how to set these environment variables permenently for a specific conda environment can be found [here](https://github.com/cmt-dtu-energy/MagTense/tree/master/python#linux-set-ld_library_path-for-specific-conda-environment-only).
    ```sh
    export LD_PRELOAD=/CONDA_ENV_PATH/lib/libmkl_core.so:/CONDA_ENV_PATH/lib/libmkl_intel_lp64.so:/CONDA_ENV_PATH/lib/libmkl_intel_thread.so:/CONDA_ENV_PATH/lib/libiomp5.so
    export LD_LIBRARY_PATH=/CONDA_ENV_PATH/lib/
    ```


## Current code development
The main features being worked on at the moment are:
- Proper code documentation
- Non-uniform grids
