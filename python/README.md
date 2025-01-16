# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the interface file `MagTense/python/src/magtense/lib/FortranToPythonIO.f90`.

## Deployment with Conda (Intel architectures)

### Requirements

- Python >= 3.12
  ```
  conda create -y -n magtense-env python==3.12 numpy matplotlib meson charset-normalizer notebook h5py tqdm
  conda activate magtense-env
  ```

- Required python packages for CUDA and MKL

  Available CUDA versions can be found here: [https://anaconda.org/nvidia/cuda](https://anaconda.org/nvidia/cuda)\
  *Note: Use `nvcc --version` or `nvidia-smi` to detect the correct CUDA version for your system.*

    ```bash
  conda install -y -c "nvidia/label/cuda-12.6.3" cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev libnvjitlink-dev
  ```

  More information about the Intel Compilers: [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

  - [ Linux ]
    ```bash
    conda install -y -c https://software.repos.intel.com/python/conda/ -c conda-forge mkl mkl-devel mkl-static "dpcpp_linux-64" intel-fortran-rt "ifx_linux-64" ncurses
    ```

  - [ Windows ] Additonal installation of Make utility

    ```bash
    conda install -y -c https://software.repos.intel.com/python/conda/ -c conda-forge mkl mkl-devel mkl-static "dpcpp_win-64" intel-fortran-rt "ifx_win-64" make
    ```

#### Compile Fortran source files
- Linux
  ```bash
  cd MagTense/python/src/magtense/lib/
  make
  ```

- Windows

  - Installation of [Visual Studio 2022](https://visualstudio.microsoft.com)

  - When `make` can not be found during Makefile execution, change to Developer PowerShell.
    This can be added to your profiles in VS Code by adding the follwing to `settings.json`:

    ```bash
    "Developer PowerShell for VS 2022": {
            "source": "PowerShell",
            "icon": "terminal-powershell",
            "args": [
              "-NoExit",
              "-ExecutionPolicy",
              "ByPass",
              "-File",
              "C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/Tools/Launch-VsDevShell.ps1"
            ]
        },
    ```
    ```bash
    cd MagTense/python/src/magtense/lib/
    make ps
    ```

  - Compilation with `nvcc` should be executed in `x64 Native Tools Command Prompt for VS 2022`.
    Otherwise, `x86` will be silently used, which results in `error: asm operand type size(8) does not match type/size implied by constraint 'r'` in `cuda_bf16.hpp`.
    Also, `f2py` needs to be run in that shell to make `ifx` compiler available for meson. 
    The x64 Command Prompt can be added to your profiles in VS Code by adding the follwing to `settings.json`:

    ```bash
    "DevCmd": {
          "path": [
              "${env:windir}\\Sysnative\\cmd.exe",
              "${env:windir}\\System32\\cmd.exe"
          ],
          "args": [
              "/d",
              "/k", 
              "C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvarsall.bat",
              "amd64"],
          "icon": "terminal-cmd"
      },
      ```

      ```bash
      cd MagTense/source/MagTenseFortranCuda/cuda
      make ps
      cd MagTense/python/src/magtense/lib/
      make 
      ```

### Installation from source (including MacOS Intel & ARM)

Create an importable Python module from Fortran source code.

For MacOS ARM architectures, currently only magnetostatics with the gfortran compiler is supported.
Navigate to folder `MagTense/python/src/magtense/lib/`, run `make`, and install the package:

```bash
cd MagTense/python/
python -m pip install -e .
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


## Distribution on [Anaconda](https://anaconda.org/cmt-dtu-energy/magtense)

### Required compilers have to be pre-installed

- [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#inpage-nav-6-undefined)
- [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

```bash
conda install -y anaconda-client conda-build

# Add nvidia channel to find CUDA and intel libraries
# conda config --show channels
conda config --env --append channels nvidia/label/cuda-12.6.3
conda config --env --append channels https://software.repos.intel.com/python/conda/
conda config --env --append channels conda-forge

# Quick fix for now
# Copy pre-compiled Python extension to MagTense/python/src/magtense/lib
# Build conda seperately for each version

# Version numbers have to be set in advance in pyproject.toml
cd MagTense/python/
python scripts/dist_conda.py
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
