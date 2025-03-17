<div align="center">
  <picture>
    <source media="(prefers-color-scheme: light)" srcset="./docs/source/static/MagTense_logo.png" height=250>
    <img alt="MagTense Logo" src="./docs/source/static/MagTenseLogo_White.svg" height=250>
  </picture>
  <br>
</div>

# MagTense

MagTense is a framework for magnetostatic and micromagnetic calculations.

## Features

- Interfaces for MATLAB and Python, with the core implemented in Fortran for speed;
- Fully analytical calculation of demagnetization tensor for cylinders, pieces of cylinders, prisms, circular pieces and tetrahedrons;
- Micromagnetic solutions of the Landau-Lifshitz equations, using the analytical demagnetization tensor described above;
- GPU-accelerated code with [CUDA®](https://developer.nvidia.com/cuda-zone) (requires NVIDIA graphics card).
- Uses Intel MKL for the micromagnetic simlations and can also utilize CVODE.
- Tested in Linux and Windows 11+ (macOS not supported at the moment).

## Installation and usage with the Python interface

Installation is recommended via `conda` package manager (requires >=**Python 3.9**).

- Installation with CUDA 12 (requires NVIDIA graphics card):
  
  ```
  conda install magtense -c cmt-dtu-energy/label/cuda-12 -c nvidia/label/cuda-12.6.3 -c https://software.repos.intel.com/python/conda/ -c conda-forge
  ```

- Installation without CUDA support:

  ```
  conda install magtense -c cmt-dtu-energy/label/cpu -c https://software.repos.intel.com/python/conda/ -c conda-forge
  ```

Examples of how to calculate magnetostatic and micromagnetic problems using the Python interface can be found [python/examples/](./python/examples/).

### Building from source

If you want to build MagTense yourself, you can find detailed instructions for building the Python interface with the Fortran core [here](/python/README.md).

## Installation and usage with the MATLAB interface

MagTense is directly useable in Matlab on Windows by downloading the MEX-files in [Releases](https://github.com/cmt-dtu-energy/MagTense/releases). Only Matlab 2020b or greater is required.

Examples of how to calculate magnetostatic and micromagnetic problems using the Matlab interface can be found in [matlab/examples](./matlab/examples).

### Building from source

If you want to compile MagTense on Windows, we provide a Visual Studio project file, [MagTense.sln](MagTense.sln), while for Linux we provide a Makefile. To build the MEX-files, we provide a Matlab function called [buildMagTenseMEX.m](matlab/buildMagTenseMEX.m) that works on both OS. You can find more information [here](/matlab/README.md)

## Further documentation

The webpage of the code is available at https://www.magtense.org.

The TechManual on the code is available at https://cmt-dtu-energy.github.io/MagTense.

## Citation
If you use this package in a publication, or simply want to refer to it, please cite the paper below:

```bibtex
@article{BJORK2021168057,
title = {MagTense: A micromagnetic framework using the analytical demagnetization tensor},
journal = {Journal of Magnetism and Magnetic Materials},
volume = {535},
pages = {168057},
year = {2021},
issn = {0304-8853},
doi = {https://doi.org/10.1016/j.jmmm.2021.168057},
url = {https://www.sciencedirect.com/science/article/pii/S0304885321003334},
author = {R. Bj{\o}rk and E. B. Poulsen and K. K. Nielsen and A. R. Insinga},
}
