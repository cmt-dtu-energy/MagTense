<div align="center">
  <img src="https://cmt-dtu-energy.github.io/MagTense/_static/MagTense_logo.PNG" height=250><br>
</div>

# MagTense version 2023

MagTense consists of both a magnetostatic and a micromagnetism calculation framework.

The magnetostatic framework can calculate the magnetic field from objects shaped as cylinders, pieces of cylinders, prisms, circular pieces and tetrahedrons. This is done using a fully analytical calculation of the demagnetization tensor. The framework is fully implemented in Fortran and has both a Matlab MEX interface and a Python interface.

The micromagnetism framework solves the Landau-Lifshitz equation. The framework is fully implemented in Fortran and has a Matlab MEX interface and a Python interface, as well as an older Matlab implementation. The micromagnetism framework utilizes the magnetostatic framework for calculating the demagnetization field.

The webpage of the code is available at https://www.magtense.org.

The TechManual on the code is available at https://cmt-dtu-energy.github.io/MagTense.


## Usage with Matlab

MagTense is directly useable in Matlab on Windows by downloading the already compiled MEX-files in [Releases](https://github.com/cmt-dtu-energy/MagTense/releases). The files are directly useable with no compilation required, although Matlab 2020b or greater is required. Examples of how to calculate magnetostatic and micromagnetic problems using the Matlab interface can be found in [matlab/examples](https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples).


### Compilation with a Visual Studio project file

If you want to compile MagTense with a Visual Studio project file for Windows, [MagTense.sln](https://github.com/cmt-dtu-energy/MagTense/blob/master/MagTense.sln), is available, as well as a Matlab function to build the MEX-files, [buildMagTenseMEX.m](https://github.com/cmt-dtu-energy/MagTense/blob/master/buildMagTenseMEX.m). MagTense utilizes Intel MKL for the micromagnetic simlations and can also utilize CUDA and CVODE.


## Usage with Python interface

Instructions on how to build and use the Python interface are listed in [python](https://github.com/cmt-dtu-energy/MagTense/tree/master/python). Installation is recommended via `conda` package manager (requires >=**Python 3.9**). Additionally, binary installers for the Python interface are available at the [Python Package Index (PyPI)](https://pypi.org/project/magtense).

- Installation with CUDA 11.x:
  
  ```
  conda install magtense -c cmt-dtu-energy/label/cuda-11 -c nvidia/label/cuda-11.8.0
  ```

- Installation with CUDA 12.x:
  
  ```
  conda install magtense -c cmt-dtu-energy/label/cuda-12 -c nvidia/label/cuda-12.2.2
  ```

- Installation without CUDA support:

  ```
  conda install magtense -c cmt-dtu-energy/label/cpu
  ```

## Current code development
The main features being worked on at the moment are:
- Proper code documentation
- Non-uniform grids
