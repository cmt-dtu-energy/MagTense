# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the [interface file](./FortranToPythonIO.f90).

## Deployment with Conda (Intel architectures)

### Create an importable Python module from Fortran source code

#### Linux

- New Conda environment from [yml-file](./.build/env-313-linux.yml)

  ```bash
  conda env create -n magtense-env -f python/.build/env-313-linux.yml
  conda activate magtense-env
  ````

- Required modules for `cvode` from sundials-7.2.1

  - Download version 7.2.1 of `cvode`:

    ```bash
    wget https://github.com/LLNL/sundials/releases/download/v7.2.1/cvode-7.2.1.tar.gz
    tar -xf cvode-7.2.1.tar.gz
    ```

  - Prepare folder structure

    ```bash
    mkdir cvode
    mv cvode-7.2.1 cvode/src
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
  make python USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0
  ```


  **Note:** In case the error `meson.build:1:0: ERROR: Executables created by fortran compiler ifx are not runnable.` shows up, most probably, some shared object dependencies from the Conda environment cannot be found.
  A temporary solution to this is adding the library path manually to `LD_LIBRARY_PATH`:

    ```bash
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
    ```

#### Windows

- Conda environment from `env-313-win.yml`

  ```bash
  conda env create -f python/.build/env-313-win.yml
  ```

  OR

  ```bash
  conda create -y -n magtense-env
  conda activate magtense-env
  conda config --env --add channels conda-forge
  conda install -y python=3.13
  python -m pip install numpy meson charset-normalizer
  conda config --env --add channels https://software.repos.intel.com/python/conda/
  conda install -y mkl mkl-devel mkl-static "dpcpp_win-64" intel-fortran-rt "ifx_win-64"
  conda config --env --add channels nvidia/label/cuda-12.8.1
  conda install -y cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev libnvjitlink-dev
  conda install -y git make
  ```

- Required modules for `cvode` from sundials-7.2.1

  - Installation of [Visual Studio 2022](https://visualstudio.microsoft.com) (Desktop development with C++ and `cmake`)

  - Installation of [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) (both C++ and Fortran)

  - Download [sundials-7.2.1](https://github.com/LLNL/sundials/releases/download/v7.2.1/cvode-7.2.1.tar.gz) and unzip it.

  Open a "Intel oneAPI command prompt for Intel 64 for Visual Studio 2022" as administrator and then do:
  ```bash
  "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -G "Ninja" -B C:/CVODE_temp -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_Fortran_COMPILER=ifx -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DENABLE_OPENMP=ON

  cmake --build C:/CVODE_temp --config RELEASE --verbose
  cmake --install C:/CVODE_temp --verbose
  ```
  where the `CVODE_temp` is a temporary directory that can be removed after installation.
  *Note: The commands below assume the `Enterprise` version of Visual Studio being installed - if you have the free community edition, simply change `Enterprise` to `Community` in the path.*

  - The installed CVODE files will be located in `"C:\Program Files (x86)\SUNDIALS"` but should be moved to a folder named `cvode` at the top-level of the MagTense repository:

    ```bash
    mkdir cvode
    xcopy "C:\Program Files (x86)\SUNDIALS\*" cvode /s /i
    ```

- Compile Fortran source files
  
  - Installation of [Visual Studio 2022](https://visualstudio.microsoft.com) / Desktop development with C++

  - Setting up customized versions of `Developer PowerShell` and `x64 Native Tools Command Prompt for VS 2022`:

    - In [VS Code](https://code.visualstudio.com), these integrated terminals can be added to your profiles by editing `settings.json`. Further, the Python extension ensures that the correct Conda environment is activated in all terminals.

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
    make magnetostatic micromagnetism USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0
    ```

  - Compilation with `nvcc` should be executed in `x64 Native Tools Command Prompt for VS 2022`.
    Otherwise, `x86` will be silently used, which results in `error: asm operand type size(8) does not match type/size implied by constraint 'r'` in `cuda_bf16.hpp`.

    ```bash
    cd source/MagTenseFortranCuda/cuda
    make
    ```

    - **Note:** In case the error `nvcc fatal   : Could not set up the environment for Microsoft Visual Studio [...]` shows up, the environment path in the active Conda environment prevents `nvcc` to work correctly. A quick fix to compile `MagTenseCudaBlas` is to initialize a `x64 Native Tools Command Prompt for VS 2022` without `conda`:

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
    make python-win USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0
    ```

    - **Note:** In case error `meson.build:1:0: ERROR: Unknown compiler(s): [['ifx']]` shows up, it should help to reinitialize your Conda environment to ensure having the correct environment path:

      ```bash
      conda deactivate
      conda activate magtense-env
      ```

### Install local editable magtense package

```bash
cp python/.build/requirements-py3-dev.txt python/requirements.txt
python3 -m pip install -e ./python
```

### Required packages at runtime

The `python/.build/` contains requirement-files, which are shipped with the respective pip-wheel.

```bash
python3 -m pip install numpy mkl intel-fortran-rt matplotlib notebook h5py tqdm importlib_resources
python3 -m pip install nvidia-cuda-runtime-cu12 nvidia-cublas-cu12 nvidia-cusparse-cu12 nvidia-nvjitlink-cu12 # only required for cuda support
```

### Latest versions of required packages

- CUDA

  Available CUDA versions can be found here: [https://anaconda.org/nvidia/cuda](https://anaconda.org/nvidia/cuda) \
  Location of corresponding [https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#pip-wheels](pip-wheels) for deployment \
  *Note: Use `nvcc --version` or `nvidia-smi` to detect the correct CUDA version for your system.*

  ```bash
  conda config --env --add channels nvidia/label/cuda-12.9.1
  conda install -y cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev libnvjitlink-dev
  ```

- Intel compilers and MKL

  More information about the Intel Compilers: [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

  ```bash
  conda config --env --add channels https://software.repos.intel.com/python/conda/
  conda install -y mkl mkl-devel mkl-static "dpcpp_linux-64" "ifx_linux-64"
  ```