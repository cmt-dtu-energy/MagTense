# Using MagTense with Matlab

Using MagTense with Matlab is as easy as downloading the latest [release](https://github.com/cmt-dtu-energy/MagTense/releases) as this contains MEX-files for both Windows and Linux. All that is required is Matlab 2023a or later and an installation of the CUDA toolkit if you want to run with CUDA. The example scripts provided with the release, also available [here](https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples), shows how to call the MEX-functions from Matlab.

# Compiling MEX-files for Matlab yourself

If you want to compile your own MEX-files, the guide below shows how to do so. Building the MEX-files is always a two-step procedure. First the Fortran objectives files must be compiled and then the Matlab MEX compiler must be used to build the MEX-files. Included in MagTense is a Matlab function to build the MEX-files, [buildMagTenseMEX.m](buildMagTenseMEX.m). MagTense utilizes Intel MKL for the micromagnetic simlations and can also utilize CUDA and CVODE. 

## Compilation on Linux 

### Using make
To build the MEX-files on Linux using make (at present without the additional CVODE libraries), one can simply use the python environment provided with MagTense to install all required compilers etc. All that is required is an installation of Matlab present on the machine.

One starts by doing a clone of the MagTense repository, followed by installing the required linux packages, which are
```
sudo apt install make
sudo apt install gcc
sudo apt install gfortran
```
where gfortran is only required to make Matlab acknowledge that a fortran MEX-compiler exists - otherwise it is not utilized.

Next Miniconda (or Anaconda) must be installed, after which the MagTense conda environment with all necessary compilers can be created using

```
conda env create -f python/.build/env-313-linux.yml
conda activate magtense-env
```

Following this compilation of the MagTense Fortran files is as simple as

```
make magnetostatic micromagnetism cuda forceintegrator USE_CUDA=1 USE_CVODE=0 USE_MATLAB=1 MATLAB_INCLUDE=path_to_matlab/extern/include MKL_ROOT=${CONDA_PREFIX}
```
Please set the `path_to_matlab` and set `USE_CUDA=0` if the machine does not have a CUDA-supported GPU.

Finally, start up Matlab (If you have installed Intel oneAPI and not used the Conda environment, you will need to do `. /opt/intel/oneapi/setvars.sh` first), and then setup the MEX-compiler using `mex -setup FORTRAN` and then run 
```
buildMagTenseMEX('USE_RELEASE', true, 'USE_CUDA', false, 'USE_CVODE', false)
```
to build the MEX-files. If you want to build without CUDA, please do `buildMagTenseMEX('USE_RELEASE', false, 'USE_CUDA', false, 'USE_CVODE', false)`


## Compilation on Windows
On Windows there are two compilation options for building the Fortran objective files - using Visual Studio or using `make` through conda. Both are described below. After either step, the  [buildMagTenseMEX.m](buildMagTenseMEX.m) function must be used to build the MEX-files.

#### Regarding the CUDA objective files
If you want to compile MagTense with CUDA, an initial step is necessary. In Windows the CUDA object files cannot be build directly be neither Visual Studio nor by `make`, but must be compile in two seperate steps. The reason for this is that Intel's Fortran Compiler and the PG Fortran compiler (which is used for doing CUDA directly in Fortran) are not compatible. One simply cannot link object files from the two compilers as there is no standard for this. The strategy then is to make the CUDA GPU kernel in C++, compile this with the nvcc compiler (that uses MSVC at the core), output to an object file and then compile a C++ wrapper with Intel C++ compiler (icx) that includes and uses the nvcc compiled file. The output here should be an `obj` file that can be called from Fortran (using the Intel Fortran compiler) via the standard way (iso_c_binding) for calling C++ functions from Fortran.

### Compilation with Visual Studio

#### CUDA step 1

To compile the CUDA objective files, you need an installation of Visual Studio. The commands below assume that you have the Enterprise version of Visual Studio - if you have the free community edition, simply change `Enterprise` to `Community` in the path.

From a standard command prompt navigate to the MagTense dir "MagTense/source/MagTenseFortranCUDA/cuda" and compile with nvcc, using a specific CUDA version and the latest version of Visual Studio: 

```bash
"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.XX\bin\nvcc.exe" -c MagTenseCudaBlas.cu -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.43.34808\bin\Hostx64\x64"
```

#### CUDA step 2

Start an Intel 64 prompt (Start Menu -> Intel oneAPI 20XX -> Intel oneAPI command prompt for Intel 64 for Visual Studio XX) and compile the C++ wrapper with icx including the cuda stuff:

```bash
icx -c MagTenseCudaBlasICLWrapper.cxx
```

#### Compiling the objective files
For compiling MagTense with Visual Studio, we provide a VS project file for Windows, [MagTense.sln](../MagTense.sln). The Visual Studio environment has configuration for Release, Debug as well as for configurations included NO_CUDA and NO_CVODE. Remember to set the correct paths to Matlab and CUDA in `properties` in each project in Visual Studio. An installation of Intel oneAPI is necessary to compile using Visual Studio.

#### Building the MEX-files
Once the objective files have been build, start up Matlab, and first setup the MEX-compiler using `mex -setup FORTRAN` and then run 
```
buildMagTenseMEX('USE_RELEASE', true, 'USE_CUDA', true, 'USE_CVODE', false)
```
to build the MEX-files. If you want to build without CUDA, please do `buildMagTenseMEX('USE_RELEASE', false, 'USE_CUDA', false, 'USE_CVODE', false)`


### Compilation with make

Compiling with make requires installation of Intels oneAPI compilers and MKL. This can be a cumbersome process, so we therefore provide a conda environment file similar to Linux. Installing this the python environment provided with MagTense to install all required compilers etc. All that is required is an installation of Matlab present on the machine.

First of all Miniconda (or Anaconda) must be installed, after which the MagTense conda environment with all necessary compilers can be created using

```
conda env create -f python/.build/environment_win.yml
conda activate magtense-env
```

#### CUDA step 1

To compile the CUDA objective files you need an installation of Visual Studio as well. The commands below assume that you have the Enterprise version of Visual Studio - if you have the free community edition, simply change `Enterprise` to `Community` in the path.

Open an Anaconda prompt and activate the MagTense environment. Then navigate to the MagTense dir "MagTense/source/MagTenseFortranCUDA/cuda" and compile with nvcc, using a specific CUDA version and the latest version of Visual Studio: 

```bash
"nvcc.exe -c MagTenseCudaBlas.cu -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.43.34808\bin\Hostx64\x64"
```

#### CUDA step 2

Open a powershell with ExecutionPolicy bypass (`pwsh -NoExit -ExecutionPolicy ByPass`) and then do

```bash
cd "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\Common7\Tools\"
.\Launch-VsDevShell.ps1
cd MagTense/source/MagTenseFortranCUDA/cuda
make wrap
```

#### Compiling the objective files
Compiling the MagTense Fortran objective files is as simple as:

```
make magnetostatic micromagnetism forceintegrator USE_CUDA=1 USE_CVODE=0 USE_MATLAB=1 MATLAB_INCLUDE=path_to_matlab/extern/include MKL_ROOT=${CONDA_PREFIX}
```
Please set the `path_to_matlab` and set `USE_CUDA=0` if the machine does not have a CUDA-supported GPU.

#### Building the MEX-files
The next step is building the MEX-files. However, if you want to use only the Conda MagTense environment and not install Intel oneAPI, Matlab cannot find the required ifx compiler to build the MEX-files, even though it is present in the Conda environment. This can be remedied by making a symbolic link to the ifx compiler in the Conda enviroment so that Matlab's "mex -setup FORTRAN" can find it. Using a Powershell with the correct `path_to_conda_environments` do:

```bash
cd path_to_conda_environments\envs\magtense-env\Library 
mkdir compiler 
cd path_to_conda_environments\envs\magtense-env\Library\compiler
mkdir 2024.2 
cd 2024.2 
cmd /c mklink /D bin path_to_conda_environments\envs\magtense-env\Library\bin
```

Now start up Matlab. To find the Intel ifx compiler in the Conda MagTense environment, we set the `ONEAPI_ROOT` environment variable and then configure the MEX-compiler using `mex -setup FORTRAN` and then build the MEX-files specifying where the intel compilers are as options:
```
setenv('ONEAPI_ROOT',"C:\Users\runneradmin\miniconda3\envs\magtense-env\Library");
mex -setup FORTRAN;
cd('MagTense/matlab');
buildMagTenseMEX('USE_CUDA',true,'USE_CVODE',false,'mkl_include','path_to_conda_environments\envs\magtense-env\opt\compiler\include\intel64','mkl_lib','path_to_conda_environments\envs\magtense-env\Library\lib','mkl_lp64','path_to_conda_environments\envs\magtense-env\Library\include\intel64\lp64');
```

## Install CVODE from sundials-4.1.0

- Requirements:

  - [cmake](https://cmake.org/)

  - [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

  - [Sundials-4.1.0](https://github.com/LLNL/sundials/releases/tag/v4.1.0)

      ```bash
      wget https://github.com/LLNL/sundials/releases/download/v4.1.0/sundials-4.1.0.tar.gz
      tar -xf sundials-4.1.0.tar.gz
      ```

    Note: It is not possible to get the cmake installer to work if sundials is unzipped in a directory with spaces.

- Rename the unpacked folder `sundials-4.1.0` to `srcdir`

    ```bash
    mv sundials-4.1.0 srcdir
    mkdir sundials-4.1.0
    mv srcdir ./sundials-4.1.0
    cd sundials-4.1.0
    mkdir builddir
    mkdir instdir
    ```

- Open the file "scrdir\config\SundialsFortran.cmake" and change the following two lines to:

    ```bash
    214     "FIND_FILE(FLIB flib.f ${FortranTest_DIR})\n"
    262     "FIND_FILE(FLIB flib.f ${FortranTest_DIR})\n"
    ```

## Build with cmake on Windows

- Open Visual Studio and open a Developer Command Prompt (Tools/Command Line/Developer Command Prompt)

    ```bash
    cd builddir
    "C:\Program Files\CMake\bin\cmake-gui.exe" ../srcdir
    ```

- Click "Configure"

    Choose Visual Studio 16 2019 (or whatever latest version you have) as generator
    Write "Intel C++ Compiler 19.0" (or whatever latest version you have) in "Optional toolset to use"

- Set the CMAKE_INSTALL_PREFIX to the full path to the instdir
- Set the EXAMPLES_INSTALL_PATH to the full path to the instdir/examples
- Enable only the following settings (if not all settings appear, select the appropriate onces and click "Configure" again):

    ```bash
    BUILD_CVODE
    BUILD_STATIC_LIBS
    BUILD_TESTING
    EXAMPLES_ENABLE_C
    EXAMPLES_ENABLE_CXX
    EXAMPLES_ENABLE_F77
    EXAMPLES_ENABLE_F90
    F2003_INTERFACE_ENABLE
    F77_INTERFACE_ENABLE
    OPENMP_ENABLE
    ```

- Click "Generate" in cmake

    Note: If any paths contain spaces or parentheses, escape the offending symbols with \ and generate again. Likely culprit: BLAS libraries.

- Copy [build_rest.bat](build_rest.bat) into builddir

    ```bash
    cp matlab/build_rest.bat .
    ```

- In Visual Studio command window run the following commands:

    ```bash
    msbuild ALL_BUILD.vcxproj -property:Configuration=Debug
    msbuild INSTALL.vcxproj -property:Configuration=Debug

    msbuild ALL_BUILD.vcxproj -property:Configuration=Release
    msbuild INSTALL.vcxproj -property:Configuration=Release
    ```

- CVODE is now installed in the instdir. Move the folder `sundials-4.1.0` to the expected location:

    ```bash
    cd ..
    mv sundials-4.1.0 "C:\Program Files (x86)\sundials-4.1.0"
    ```

Note: The last four instructions *could* be put into the .bat file, but it's easier to isolate warnings / errors if they aren't.
Building ALL_BUILD causes two linker warnings because SUNDIALS gave compile flags to the linker as well as the compiler.
Building INSTALL can cause one (strange) warning, MSB8065, which is apparently a [bug](https://gitlab.kitware.com/cmake/cmake/issues/19737) of no consequence.

## Build with cmake on Linux

```bash
cd builddir
cmake --verbose -DCMAKE_INSTALL_PREFIX=/path/to/sundials-4.1.0/instdir -DEXAMPLES_INSTALL_PATH=/path/to/sundials-4.1.0/instdir/examples -DBUILD_ARKODE=OFF -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_CVODE=ON -DBUILD_STATIC_LIBS=ON -DBUILD_TESTING=ON -DCMAKE_Fortran_COMPILER=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/ifort -DEXAMPLES_ENABLE_C=ON -DEXAMPLES_ENABLE_CXX=ON -DEXAMPLES_ENABLE_F77=ON -DEXAMPLES_ENABLE_F90=ON -DF2003_INTERFACE_ENABLE=ON -DF77_INTERFACE_ENABLE=ON -DOPENMP_ENABLE=ON ../srcdir
```
