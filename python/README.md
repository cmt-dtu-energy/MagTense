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

**Note: Compiling with FMM3D backend**

To compile with FMM3D `$(MagTense-Folder)/external/FMM3D/local` must be in the `LD_Library_PATH` (replace MagTense-Folder with the actual path to the MagTense repo).
After this modify the Makefile to replace `USE_FMM3D=0` with `USE_FMM3D=1`.

#### Windows

Installation on Windows is a but more contrived because of the way Windows handles environments and privileges, so the user has to do more steps manually. Fortunately, the user who wants to customize and develop MagTense will likely handle the [pre-requisites](#pre-requisites) only once, to set up the development environment, and then the build process when the Fortran or Python parts are modified should be straightforward.

##### Pre-requisites (setting up the development environment)

First, clone this repo.

Then, install the prerequisites:

- [The Anaconda distribution](https://www.anaconda.com/download)
- [Visual Studio 2022](https://visualstudio.microsoft.com/vs/older-downloads/#visual-studio-2022-and-other-products) - select the option for desktop development with C++ and `cmake`
- [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) (both C++ and Fortran)
- [sundials-7.4.0](https://github.com/LLNL/sundials/releases/download/v7.4.0/cvode-7.4.0.tar.gz) - unzip to a folder of your choice

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

After that, you can close the Intel oneAPI command prompt.

##### Building the Fortran core and Python interface

On modern Windows distributions, the above development environment setup should create a new profile on the Windows Terminal app, called "Anaconda Powershell Prompt", with Anaconda already configured. *All commands below should be issued from this terminal, and in the MagTense directory (where you cloned the repo)*.

```shell
conda env create -f python/.build/env-314-win.yml
```

Then, activate the environment with `conda activate magtense-env`.

