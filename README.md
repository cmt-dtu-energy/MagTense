<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cmt-dtu-energy/MagTense/refs/heads/master/docs/source/static/MagTenseLogo_White.svg" height=250>
    <img alt="MagTense Logo" src="https://raw.githubusercontent.com/cmt-dtu-energy/MagTense/refs/heads/master/docs/source/static/MagTense_logo.png" height=250>
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
- Uses Intel MKL for the micromagnetic simlations and can also utilize [CVODE](https://computing.llnl.gov/projects/sundials/cvode).
- Tested in Linux and Windows 11+ (macOS not supported at the moment).

## Installation and usage with the Python interface

Installation is recommended via `pip` (requires >=**Python 3.12**):

```
pip install magtense
```

Examples of how to calculate magnetostatic and micromagnetic problems using the Python interface can be found [python/examples/](./python/examples/).

## Installation and usage with the Matlab interface

MagTense is directly useable in Matlab on Windows by downloading the MEX-files in [Releases](https://github.com/cmt-dtu-energy/MagTense/releases). Only Matlab 2023a or greater is required.

Examples of how to calculate magnetostatic and micromagnetic problems using the Matlab interface can be found in [matlab/examples](./matlab/examples).

## Building from source

If you want to build MagTense yourself this is certainly also an option.

If you want to compile the MagTense core on Linux we provide a Makefile, which also works on Windows, where we also provide a Visual Studio project file, [MagTense.sln](MagTense.sln).

For the higher-level interfaces, the instructions depend on which one you want to build.

### Python interface

#### Linux

Provided you already have `make` installed, you can simply run

```shell
make python-interface
```

This will:

- Download [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) if you don't have it already (by default it assumes the installation lives at `~/miniconda3`);
- Create the conda environment `magtense-env` with all the dependencies for building the Python interface;
- Build the MagTense Fortran core and the CVODE library;
- Build the Python interface and install it in the `magtense-env` conda environment

There are also Make targets `rm-env` and `rm-conda` to clean up your installation. The Make rules should be smart enough to not do unnecessary work; for instance, if you modify some part of the Python interface, `make python` will re-build the wheel and re-install it on the environment, without re-downloading Miniconda. In case of errors, the simplest first step probably is to run `make rm-env` and `make rm-conda` to start from a clean slate.

To test if the interface was built correctly, run `make test`.

Note that all automated commands are run by subshells managed by Make. To actually use conda and explore the Python files, you have to first initialize conda in your current shell:

```shell
$HOME/miniconda3/bin/conda init --all
```

Then, *create a new shell*, and activate the `magtense-env` conda environment:

```shell
conda activate magtense-env
```

As a starting point, you can run the example scripts in [python/examples/](./python/examples/).

#### Windows

Installation on Windows is a but more contrived because of the way Windows handles environments, so the user has to do more steps manually.

##### Pre-requisites

First, clone this repo.

Then, install the prerequisites:

- [The Anaconda distribution](https://www.anaconda.com/download)
- [Visual Studio 2022](https://visualstudio.microsoft.com/vs/older-downloads/#visual-studio-2022-and-other-products) - select the option for desktop development with C++ and `cmake`
- [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) (both C++ and Fortran)
- [sundials-7.4.0](https://github.com/LLNL/sundials/releases/download/v7.4.0/cvode-7.4.0.tar.gz) - unzip to a folder of your choice

##### Building the cvode dependency



##### Building the Fortran core

##### Building the Python interface

On modern windows distributions, this should create a new profile on the Windows Terminal app, with Anaconda set up. *All commands below should be issued from this terminal, and in the MagTense directory (where you cloned the repo)*.

```shell
conda env create -f python/.build/env-314-win.yml
```

Then, activate the environment with `conda activate magtense-env`.


### MATLAB interface

For Matlab MEX-files, we provide a Matlab function called [buildMagTenseMEX.m](matlab/buildMagTenseMEX.m) that works on both OS. You can find more information [here](/matlab/README.md).

## Further documentation

The webpage of the code is available at <https://www.magtense.org>.

The TechManual on the code is available at <https://cmt-dtu-energy.github.io/MagTense>.

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
