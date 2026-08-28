# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the [interface file](./FortranToPythonIO.f90).

## Linux

On a Debian/Ubuntu system, first install the build prerequisites:

```shell
sudo apt-get update
sudo apt-get install -y git curl build-essential
```

Then you can simply run

```shell
make python-interface [PY_VERSION=314(default) | 313 | 312] [USE_CUDA=1(default) | 0]
```

This will:

- Download [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) if you don't have one already. An existing installation is detected automatically, in this order: a `CONDA_DIR` (or `CONDA_BIN`) you pass yourself, the `$CONDA_EXE` of an activated conda shell, `conda info --base` from your `PATH`, and finally `~/miniconda3`. Only when none of those turns up a conda is Miniconda downloaded, into `~/miniconda3`. To point the build at a specific installation, pass it explicitly:

  ```shell
  make python-interface CONDA_DIR=/path/to/your/conda
  ```

  You can check what was picked up with `make -s print-CONDA_DIR`, `make -s print-CONDA_BIN` and `make -s print-PYTHON` (or all of it at once with `make info`);

  The supported distributions are **Miniconda and Anaconda**: the conda front-end is assumed to be at `$(CONDA_DIR)/bin/conda`. Other front-ends (mamba, micromamba) are untested, but you can point the build at one with `CONDA_BIN=/path/to/mamba`.

  Note that `CONDA_DIR` is not remembered between invocations - every target has to be given the same value, so a later `make pytest` from a plain shell needs `make pytest CONDA_DIR=/path/to/your/conda` too. If you always build against the same installation, export `CONDA_DIR` from your `~/.bashrc` instead of passing it each time.
- Create the conda environment `magtense-env` with all the dependencies for building the Python interface, using the specified Python vesrsion
- Build the MagTense Fortran core and the CVODE library, using CUDA by default (if you have an NVIDIA GPU and the CUDA toolkit installed) - set `USE_CUDA=0` to disable CUDA support
- Build the Python interface and install it in the `magtense-env` conda environment

There are also Make targets `rm-env` and `rm-conda` to clean up your installation. `rm-conda` only removes a Miniconda that this Makefile installed itself; it refuses to touch a conda that was already on your machine. The Make rules should be smart enough to not do unnecessary work; for instance, if you modify some part of the Python interface, `make python` will re-build the wheel and re-install it on the environment, without re-downloading Miniconda. In case of errors, the simplest first step probably is to run `make rm-env` and `make rm-conda` to start from a clean slate.

To test if the interface was built correctly, run `make pytest`.

Note that all automated commands are run by subshells managed by Make. To actually use conda and explore the Python files, you have to first initialize conda in your current shell:

```shell
$(make -s print-CONDA_BIN) init --all
```

Then, *create a new shell*, and activate the `magtense-env` conda environment:

```shell
conda activate magtense-env
```

As a starting point, you can run the example scripts in [python/examples/](./python/examples/).

**Note: Compiling with FMM3D backend**

To compile with FMM3D `$(MagTense-Folder)/external/FMM3D/local` must be in the `LD_Library_PATH` (replace MagTense-Folder with the actual path to the MagTense repo).
After this modify the Makefile to replace `USE_FMM3D=0` with `USE_FMM3D=1`.

## Windows

Installation on Windows is a but more contrived because of the way Windows handles environments and privileges, so the user has to do more steps manually. Fortunately, the user who wants to customize and develop MagTense will likely handle the [pre-requisites](#pre-requisites-setting-up-the-development-environment) only once, to set up the development environment, and then the build process when the Fortran or Python parts are modified should be straightforward.

##### Pre-requisites (setting up the development environment)

First, clone this repo.

Then, install the prerequisites:

- [The Anaconda distribution](https://www.anaconda.com/download)
- [Visual Studio 2022](https://visualstudio.microsoft.com/vs/older-downloads/#visual-studio-2022-and-other-products) - select the option for desktop development with C++ and `cmake`
- [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) (both C++ and Fortran)
- [sundials-7.4.0](https://github.com/LLNL/sundials/releases/download/v7.4.0/cvode-7.4.0.tar.gz) - unzip to a folder of your choice
- *Only if you want the FMM3D backend*: [FMM3D](https://github.com/Ximtecs/FMM3D), unzipped into `external\FMM3D`

With the above install, there's a new application in the start menu called "Intel oneAPI command prompt for Intel 64 for Visual Studio 2022" - this is a terminal with all the necessary environment variables set up to use the Intel compilers and tools. Find it and open it as administrator (right-click -> "Run as administrator"). Then, navigate to the folder where you unzipped the `cvode-7.4.0` source code and build it with the following commands:

```bash
cmake -G "Ninja" -B C:/CVODE_temp -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_Fortran_COMPILER=ifx -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DENABLE_OPENMP=ON

cmake --build C:/CVODE_temp --config RELEASE --verbose
cmake --install C:/CVODE_temp --verbose
```
where the `CVODE_temp` is a temporary directory that can be removed after installation.

(**Tip**: opening the app as described above will open a rudimentary terminal with a sub-optimal UI. In the new Windows Terminal app, you can create a new profile that runs the above command to open the "Intel oneAPI command prompt for Intel 64 for Visual Studio 2022" in a more user-friendly terminal. See [this guide](https://docs.microsoft.com/en-us/windows/terminal/customize-settings/profile-settings-guide) on how to do this. The *command* to run in this new profile is `C:\Windows\System32\cmd.exe /E:ON /K ""C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022"`, and don't forget to check the "Run this profile as administrator" option in the profile settings.)

The installed CVODE files will be located in `"C:\Program Files (x86)\SUNDIALS"` but should be moved to a folder named `cvode` at the top-level of the MagTense repository:

```bash
mkdir cvode
xcopy "C:\Program Files (x86)\SUNDIALS\*" cvode /s /i
```

This is the default location the build looks for, so `CVODE_ROOT` does not have to be passed anywhere below. After that, you can close the Intel oneAPI command prompt.

##### Optional: the CUDA plugin

The CUDA plugin is built with `nvcc` and therefore needs the MSVC toolchain rather than the conda environment. Open an `x64 Native Tools Command Prompt for VS 2022` and run:

```bash
cd source/MagTenseFortranCuda/cuda
make
```

Skip this if you intend to build with `USE_CUDA=0`.

##### Building the Fortran core and Python interface

On modern Windows distributions, the above development environment setup should create a new profile on the Windows Terminal app, called "Anaconda Powershell Prompt", with Anaconda already configured. *All commands below should be issued from this terminal, and in the MagTense directory (where you cloned the repo)*.

```shell
conda env create -f python/.build/env-314-win.yml
```

Then, activate the environment with `conda activate magtense-env`.

With the environment active, build the Fortran core and then the Python module:

```shell
make auxmt magnetostatic micromagnetism USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0
make python-win USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0
```

Pass `USE_CUDA=0` if you skipped the CUDA plugin, and add `USE_FMM3D=1` (after `make fmm3d USE_FMM3D=1`) if you unzipped FMM3D. The same flags have to be given to both commands.

Finally, install the compiled package into the environment so that simulations can be run:

```shell
cp python/.build/requirements-py3-dev.txt python/requirements.txt
python -m pip install -e ./python
```

As a starting point, you can run the example scripts in [python/examples/](./examples/).
