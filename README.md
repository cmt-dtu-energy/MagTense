<div align="center">
  <picture>
    <source media="(prefers-color-scheme: light)" srcset="./docs/source/static/MagTense_logo.png" height=250>
    <img alt="MagTense Logo" src="./docs/source/static/MagTenseLogo_White.svg" height=250>
  </picture>
  <br>
</div>

# MagTense

MagTense consists of both a magnetostatic and a micromagnetism calculation framework.

The magnetostatic framework can calculate the magnetic field from objects shaped as cylinders, pieces of cylinders, prisms, circular pieces and tetrahedrons. This is done using a fully analytical calculation of the demagnetization tensor. The framework is fully implemented in Fortran and has both a Matlab MEX interface and a Python interface.

The micromagnetism framework solves the Landau-Lifshitz equation. The framework is fully implemented in Fortran and has a Matlab MEX interface and a Python interface. The micromagnetism framework utilizes the magnetostatic framework for calculating the demagnetization field.

The webpage of the code is available at https://www.magtense.org.

The TechManual on the code is available at https://cmt-dtu-energy.github.io/MagTense.


## Usage with Matlab

MagTense is directly useable in Matlab on Windows by downloading the already compiled MEX-files in [Releases](https://github.com/cmt-dtu-energy/MagTense/releases). The files are directly useable with no compilation required, although Matlab 2020b or greater is required. Examples of how to calculate magnetostatic and micromagnetic problems using the Matlab interface can be found in [matlab/examples](matlab/examples).


### Compilation with a Visual Studio project file

If you want to compile MagTense with a Visual Studio project file for Windows, [MagTense.sln](MagTense.sln), is available, as well as a Matlab function to build the MEX-files, [buildMagTenseMEX.m](matlab/buildMagTenseMEX.m). MagTense utilizes Intel MKL for the micromagnetic simlations and can also utilize CUDA and CVODE.


## Usage with Python interface

Instructions on how to build and use the Python interface are listed in [python](python). Installation is recommended via `conda` package manager (requires >=**Python 3.9**). Additionally, binary installers for the Python interface are available at the [Python Package Index (PyPI)](https://pypi.org/project/magtense).

- Installation with CUDA 12:
  
  ```
  conda install magtense -c cmt-dtu-energy/label/cuda-12 -c nvidia/label/cuda-12.6.3 -c https://software.repos.intel.com/python/conda/ -c conda-forge
  ```

- Installation without CUDA support:

  ```
  conda install magtense -c cmt-dtu-energy/label/cpu -c https://software.repos.intel.com/python/conda/ -c conda-forge
  ```

