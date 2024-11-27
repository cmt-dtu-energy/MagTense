# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the interface file `MagTense/python/src/magtense/lib/FortranToPythonIO.f90`.

## Deployment with Conda (Intel architectures)

### Requirements

- Python >= 3.9

- Required python packages

  Find available `${CUDA_LABEL}` [here](https://anaconda.org/nvidia/cuda).
  HINT: Use `nvcc --version` or `nvidia-smi` to detect the correct CUDA version for your system.

  More information about the Intel Compilers: [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)
  HINT: Set "`${OS}` to `win` or `linux`, respectively.

  ```bash
  conda install -y numpy matplotlib
  conda install -y -c "nvidia/label/cuda-${CUDA_LABEL}" cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev libnvjitlink-dev
  conda install -y -c https://software.repos.intel.com/python/conda/ -c conda-forge mkl mkl-static "dpcpp_${OS}-64" intel-fortran-rt "ifx_${OS}-64"
  conda install -y meson charset-normalizer
  ```

  - [ Windows / MacOS ] Installation of Make utility

    ```bash
    conda install -y -c conda-forge make
    ```

- Additional python packages to run scripts

    ```bash
    conda install -y notebook h5py tqdm
    ```

### Installation from source (including MacOS Intel & ARM)

Create an importable Python module from Fortran source code.

For MacOS ARM architectures, currently only magnetostatics with the gfortran compiler is supported.
Navigate to folder `MagTense/python/src/magtense/lib/`, run `make`, and install the package:

```bash
cd MagTense/python/src/magtense/lib/
make
cd MagTense/python/
python -m pip install -e .
```

## Read-in customized M-H-curve

This feature is currently only supported for soft magnetic tiles ([type=2](magtense/magtense.py#L49)).

In  [iterate_magnetization()](magtense/magtense.py#L611), an arbitrary number of state functions (M-H-curves) can be defined:

```python
mu_r = 100
datapath = f'./magtense/mat/Fe_mur_{mu_r}_Ms_2_1.csv'

 ...

data_statefcn = numpy.genfromtxt(datapath, delimiter=';')
n_statefcn = 1
```

[Here](magtense/mat), three sample M-H-curves for Fe with different relative permeabilities and a saturation magnetization of 2.1 T are stored as CSV-files. The data format is as follows:

```csv
0; Temp0; Temp1; ...
H0-field; M0@Temp0; M0@Temp1;...
H1-field; M1@Temp0; M1@Temp1;...
.
.
H100-field; M100@Temp0; M100@Temp1; ...
.
```

With only one state function given, the same M-H-curve applies to all tiles of type 2.

When the soft tiles differ in their M-H-curves, multiple state function can be combined. In order to match a specific M-H-curve with the corresponding tile, the variable [stfcn_index](magtense/magtense.py#L54) can be set.

## Distribution on [Anaconda](https://anaconda.org/cmt-dtu-energy/magtense)

### Required compilers have to be pre-installed

- [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#inpage-nav-6-undefined)
- [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

```bash
conda install -y anaconda-client conda-build

# Add nvidia channel to find CUDA and intel libraries
# conda config --show channels
conda config --env --append channels nvidia/label/cuda-12.2.2
conda config --env --append channels intel
# On Windows
# conda config --env --append channels conda-forge

# Quick fix for error of nvcc during build on Windows
# conda install -y -c nvidia/label/cuda-12.2.2 cuda-nvcc
# cd MagTense/source/MagTenseFortranCuda/cuda/
# nvcc -c MagTenseCudaBlas.cu -o MagTenseCudaBlas.o

# Version numbers have to be set in advance in pyproject.toml
cd MagTense/python/
python scripts/dist_conda.py
```

## Distribution on [PyPI](https://pypi.org/project/magtense/)

Libraries have to be pre-build for now, and should be located in `MagTense/python/compiled_libs`.

```bash
# Required python packages for distribution
python -m pip install build
conda install -y twine

cd MagTense/python/
python scripts/dist_pypi.py

# Upload to pypi.org
# twine upload --repository testpypi dist/*
twine upload dist/*
```
