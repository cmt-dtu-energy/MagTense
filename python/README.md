# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the interface file `MagTense/python/magtense/lib/FortranToPythonIO.f90`.

## Deployment with Conda (Intel architectures)

### Requirements

- Python >= 3.9

- [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#inpage-nav-6-undefined) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

  - [Linux] Prepare your terminal, so that ifort compiler can be found:

    ```bash
    . /opt/intel/oneapi/setvars.sh
    ```

  - [Windows] Activation via [VS Code extension for Intel® oneAPI Toolkits](https://github.com/intel/vscode-oneapi-environment-configurator)

- Required python packages

  Find available `${CUDA_LABEL}` [here](https://anaconda.org/nvidia/cuda).
  HINT: Use `nvcc --version` or `nvidia-smi` to detect the correct CUDA version for your system.

  ```bash
  conda install -y numpy matplotlib
  conda install -y -c "nvidia/label/cuda-${CUDA_LABEL}" cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev libnvjitlink
  conda install -y -c intel mkl-static
  ```

  - [ Linux ]

    ```bash
    conda install -y -c intel intel-fortran-rt
    ```

  - [ Windows / MacOS ] Installation of Make utility

    ```bash
    conda install -y -c conda-forge make
    ```

- Additional python packages to run data creation scripts

    ```bash
    conda install -y h5py tqdm
    ```

### Installation from source (including MacOS Intel & ARM)

Create an importable Python module from Fortran source code.

For MacOS ARM architectures, currently only magnetostatics with the gfortran compiler is supported.
Navigate to folder `MagTense/python/magtense/lib/`, run `make`, and install the package:

```bash
cd MagTense/python/magtense/lib/
make
cd MagTense/python/
pip install -e .
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

## Distribution on [PyPI](https://pypi.org/project/magtense/) and [Anaconda](https://anaconda.org/cmt-dtu-energy/magtense)

Libraries have to be pre-build for now, and should be located in `MagTense/python/magtense/compiled_libs`.

```bash
# Required conda packages for distribution
conda install -y twine anaconda-client conda-build

cd MagTense/python/

# Create wheels
python scripts/dist_pypi.py

# Upload to pypi.org
# twine upload --repository testpypi dist/*
twine upload dist/*

# Build tarball to distribute on anaconda.org
# Add nvidia channel to find CUDA libraries
conda config --show channels
conda config --env --append channels nvidia/label/cuda-${CUDA_LABEL}

cd MagTense/.conda-build
conda-build .

# Upload to anadonda
anaconda upload --user cmt-dtu-energy ${CONDA_PREFIX}/conda-bld/${ARCH}/magtense-${MT_VERSION}-py${PY}_cuda${CUDA_VERSION}.tar.bz2
```
