# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the interface file `MagTense/python/src/magtense/lib/FortranToPythonIO.f90`.

## Deployment with Conda (Intel architectures)

For MacOS ARM architectures, currently only magnetostatics with the `gfortran` compiler is supported.

### Create an importable Python module from Fortran source code

#### Linux

- New conda environment with Python >= 3.12
  ```bash
  conda create -y -n magtense-env 
  conda activate magtense-env
  conda config --env --add channels conda-forge
  conda config --env --add channels nvidia/label/cuda-12.6.3
  conda config --env --add channels https://software.repos.intel.com/python/conda/
  conda install -y python
  conda install -y numpy matplotlib meson charset-normalizer ncurses git notebook h5py tqdm
  ```

- Required python packages for CUDA and MKL

  Available CUDA versions can be found here: [https://anaconda.org/nvidia/cuda](https://anaconda.org/nvidia/cuda)\
  *Note: Use `nvcc --version` or `nvidia-smi` to detect the correct CUDA version for your system.*

  ```bash
  conda install -y cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev libnvjitlink-dev
  ```

  More information about the Intel Compilers: [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

  ```bash
  conda install -y mkl mkl-devel mkl-static "dpcpp_linux-64" intel-fortran-rt "ifx_linux-64"
  ```

- Compile Fortran source files
  ```bash
  cd python/src/magtense/lib/
  make
  ```

#### Windows

- Conda environment from `environment.yml`

  ```bash
  conda env create -f python/environment.yml
  ```

- Compile Fortran source files
  
  - Installation of [Visual Studio 2022](https://visualstudio.microsoft.com) / Desktop development with C++

  - Setting up customized versions of `Developer PowerShell` and `x64 Native Tools Command Prompt for VS 2022`:

    - In [VS Code](https://code.visualstudio.com), these integrated terminals can be added to your profiles by editing `settings.json`. Further, the Python extension ensures that the correct `conda` environment is activated in all terminals.
      ```bash
      "terminal.integrated.profiles.windows": {
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
          "DevCmdx64": {
              "path": [
                  "${env:windir}\\Sysnative\\cmd.exe",
                  "${env:windir}\\System32\\cmd.exe"
              ],
              "args": [
                  "/k",
                  "C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/Tools/VsDevCmd.bat",
                  "-startdir=None",
                  "-arch=amd64",
                  "-host_arch=x64"
              ],
              "icon": "terminal-cmd"
          },
      },
      ```

    - Alternatively, customized profiles can be set up in the `Windows Terminal` app with the following guide: https://learn.microsoft.com/en-us/windows/terminal/install#settings-json-file

      Further, `conda` has to be initialized in these terminals to have access to the `make` executable and other necessary packages.

  - Open a `Developer PowerShell` and run:

    ```bash
	conda activate magtense-env
    cd python/src/magtense/lib/
    make ps
    ```

  - Compilation with `nvcc` should be executed in `x64 Native Tools Command Prompt for VS 2022`.
    Otherwise, `x86` will be silently used, which results in `error: asm operand type size(8) does not match type/size implied by constraint 'r'` in `cuda_bf16.hpp`.

    ```bash
    cd source/MagTenseFortranCuda/cuda
    make
    ```

    - **Note:** In case the error `nvcc fatal   : Could not set up the environment for Microsoft Visual Studio [...]` shows up, the environment path in the active conda environment prevents `nvcc` to work correctly. A quick fix to compile `MagTenseCudaBlas` is to initialize a `x64 Native Tools Command Prompt for VS 2022` without `conda`:

      ```bash
      "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\bin\nvcc.exe" -c MagTenseCudaBlas.cu -o MagTenseCudaBlas.o
      ```

      And then only compile `MagTenseCudaBlasICLWrapper.cxx` before creating `libCuda` in the activated environment:

      ```bash
      cd source/MagTenseFortranCuda/cuda
      make wrap
      ```

  - Linking and wrapping libraries with `f2py` needs to be run in `x64 Native Tools Command Prompt for VS 2022` to make `ifx` compiler available for `meson`:

    ```bash
	conda activate magtense-env
    cd python/src/magtense/lib/
    make cmdx64
    ```

    - **Note:** In case error `meson.build:1:0: ERROR: Unknown compiler(s): [['ifx']]` shows up, it should help to reinitialize your conda environment to ensure having the correct environment path:

      ```bash
      conda deactivate
      conda activate magtense-env
      ```

### Install local editable magtense package

```bash
cd python/
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
