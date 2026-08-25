# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the [interface file](./FortranToPythonIO.f90).

## Linux

- New Conda environment from [env-313-linux.yml](./.build/env-313-linux.yml)

  ```bash
  conda env create -n magtense-env -f python/.build/env-313-linux.yml
  conda activate magtense-env
  ```

  OR
  
  Make your own enviroment with the required python packages for CUDA, Intel compilers (`ifx` and `icx`), `mkl` and `cmake`:
  - Available [CUDA versions](https://anaconda.org/nvidia/cuda) and location of corresponding [pip-wheels](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#pip-wheels) (deployment)
  - More information about [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

  ```bash
  conda create -y -n magtense-env && conda activate magtense-env
  conda config --env --add channels conda-forge
  conda install -y python=3.13
  python3 -m pip install numpy meson ninja charset-normalizer build
  conda config --env --add channels nvidia/label/cuda-12.9.1
  conda install -y cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev libnvjitlink-dev
  conda config --env --add channels https://software.repos.intel.com/python/conda/
  conda install -y mkl mkl-devel mkl-static "dpcpp_linux-64" "ifx_linux-64"
  conda install -y cmake
  ```

- Required modules for `cvode` from sundials-7.4.0
    Prerequisites, e.g. `ifx`, `icx` and `cmake`, have already been installed in the previous steps.

  - Download version 7.4.0 of `cvode`:

    ```bash
    wget https://github.com/LLNL/sundials/releases/download/v7.4.0/cvode-7.4.0.tar.gz
    tar -xf cvode-7.4.0.tar.gz
    ```

  - Prepare folder structure

    ```bash
    mkdir cvode
    mv cvode-7.4.0 cvode/src
    ```

  - Run `cmake` and `make`for installation

    ```bash
    cmake \
    -B cvode/build \
    -S cvode/src \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_ARKODE=OFF \
    -D BUILD_CVODE=ON \
    -D BUILD_CVODES=OFF \
    -D BUILD_IDA=OFF \
    -D BUILD_IDAS=OFF \
    -D BUILD_KINSOL=OFF \
    -D BUILD_SHARED_LIBS=OFF \
    -D BUILD_STATIC_LIBS=ON \
    -D CMAKE_INSTALL_PREFIX=cvode \
    -D EXAMPLES_INSTALL_PATH=cvode/examples \
    -D CMAKE_C_COMPILER=$(which icx) \
    -D CMAKE_Fortran_COMPILER=$(which ifx) \
    -D BUILD_FORTRAN_MODULE_INTERFACE=ON \
    -D ENABLE_OPENMP=ON
    cmake --build cvode/build --config Release --verbose
    cmake --install cvode/build --verbose
    ```

- Compile Fortran source files

  ```bash
  LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" make python USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=0
  ```

- Compiling with FMM3D backend

To compile with FMM3D `$(MagTense-Folder)/external/FMM3D/local` must be in the `LD_Library_PATH` (replace MagTense-Folder with the actual path to the MagTense destiation).
After this compile as above but with `USE_FMM3D=1`

## Windows

- Requiments

  - Installation of [Visual Studio 2022](https://visualstudio.microsoft.com) (Desktop development with C++ and `cmake`)

  - Installation of [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) (both C++ and Fortran)

  - Installation of Miniconda or Anaconda

  - Download of MagTense and [FMM3D](https://github.com/Ximtecs/FMM3D). Unzip the latter in `external\FMM3D`


- Conda environment `magtense-env` is created from `env-313-win.yml` in a `Powershell` as

  ```bash
  conda env create -f python/.build/env-313-win.yml
  conda activate magtense-env
  ```
  All commands below, except for CVODE compilation, takes place in the `magtense-env` environment

- *Optional* CUDA plugin installation.
    Open an `x64 Native Tools Command Prompt for VS 2022` and do:

    ```bash
    cd source/MagTenseFortranCuda/cuda
    make
    ```
	
- *Optional* CVODE plugin installation

  - Download [sundials-7.4.0](https://github.com/LLNL/sundials/releases/download/v7.4.0/cvode-7.4.0.tar.gz) and unzip it.
  
  - Install [choco](https://chocolatey.org/install) and then Ninja as `choco install ninja`
  
  - Open a `Intel oneAPI command prompt for Intel 64 for Visual Studio 2022` as administrator and do:
  
  ```bash
  cd cvode-7.4.0 
  mkdir install
  "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -G "Ninja" -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_INSTALL_PREFIX=install -DEXAMPLES_INSTALL_PATH=install/examples -DCMAKE_C_FLAGS=-Wno-deprecated-declarations -DCMAKE_C_COMPILER=icx-cl -DCMAKE_CXX_COMPILER=icx-cl -DCMAKE_Fortran_COMPILER=ifx -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DENABLE_OPENMP=ON
  "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" --build build --config Release --verbose
  "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" --install build --verbose
  ```
  
- Compile Fortran source files
  
  - Open a `x64 Native Tools Command Prompt for VS 2022` and run:

    ```bash
    make fmm3d USE_FMM3D=1
    make auxmt magnetostatic micromagnetism USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=1 CVODE_ROOT="cvode-7.4.0/install"
	make python-win USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=1 CVODE_ROOT="cvode-7.4.0/install"
    ```

- Install local editable magtense package
  To install the compiled MagTense package locally, so that simulations can be run, do

  ```bash
  cp python/.build/requirements-py3-dev.txt python/requirements.txt
  python -m pip install -e ./python
  ```
