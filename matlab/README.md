# Using MagTense with Matlab

Using MagTense with Matlab is as easy as downloading the latest [release](https://github.com/cmt-dtu-energy/MagTense/releases) as this contains MEX-files for both Windows and Linux. All that is required is Matlab 2023a or later and an installation of the CUDA toolkit if you want to run with CUDA. The example scripts provided with the release, also available [here](https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/examples), shows how to call the MEX-functions from Matlab.

# Compiling MEX-files for Matlab

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

- First follow the [python guide](https://github.com/cmt-dtu-energy/MagTense/tree/master/python#Windows) to install the Conda environment, CVODE and CUDA.


### Compilation with make

- Compile Fortran source files
  
  - Open a `x64 Native Tools Command Prompt for VS 2022` and run:

    ```bash
    make auxmt magnetostatic micromagnetism forceintegrator USE_CUDA=1 USE_CVODE=1 USE_MATLAB=1 USE_FMM3D=1 CVODE_ROOT="sundials-7.4.0/install" MATLAB_INCLUDE="C:/Program Files/MATLAB/R2024b/extern/include"
    ```

- Building the MEX-files
  -Next step is building the MEX-files. Start up Matlab and then build the MEX-files specifying where the intel compilers are as options:
  ```
  cd('MagTense/matlab');
  buildMagTenseMEX('USE_CUDA',true,'USE_CVODE',true,'USE_FMM',true,'VS_STUDIO',false,'USE_RELEASE',true)
  ```
  
### Compilation with Visual Studio

For compiling MagTense with Visual Studio, we provide a VS project file for Windows, [MagTense.sln](../MagTense.sln). The Visual Studio environment has configuration for Release, Debug as well as for configurations included NO_CUDA, NO_CVODE and NO_FMM. Remember to set the correct paths to Matlab and CUDA in `properties` in each project in Visual Studio. An installation of Intel oneAPI is necessary to compile using Visual Studio.

#### Building the MEX-files
Once the objective files have been build, start up Matlab, and first setup the MEX-compiler using `mex -setup FORTRAN` and then run 
```
cd('MagTense/matlab');
buildMagTenseMEX('USE_CUDA',true,'USE_CVODE',true,'USE_FMM',true,'VS_STUDIO',true,'USE_RELEASE',true)
```
to build the MEX-files.