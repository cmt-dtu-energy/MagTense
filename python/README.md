# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the [interface file](./FortranToPythonIO.f90).

## Deployment with Conda (Intel architectures)

### Create an importable Python module from Fortran source code

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

