# Usage with Matlab

## Compilation with Visual Studio

If you want to compile MagTense with a Visual Studio project file for Windows, [MagTense.sln](../MagTense.sln), is available, as well as a Matlab function to build the MEX-files, [buildMagTenseMEX.m](buildMagTenseMEX.m). MagTense utilizes Intel MKL for the micromagnetic simlations and can also utilize CUDA and CVODE. The Visual Studio environment has configuration for Release, Debug as well as for configurations included NO_CUDA and NO_CVODE.

## Compilation with make

### Requirements

- [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#inpage-nav-6-undefined) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

- **Intel® MKL** [directly](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) or via [conda](https://anaconda.org/intel/mkl-static) | Expected location:

  [ Windows ] `C:/Program Files (x86)/Intel/oneAPI/mkl/latest`

  [ Linux + oneapi ] `/opt/intel/oneapi/mkl/latest`

  [ Linux + conda ] `$CONDA_PREFIX/lib/intel64`

- [ Linux + conda ] [Runtime for Intel® Fortran Compiler](https://anaconda.org/intel/intel-fortran-rt)

- [ Windows / MacOS ] **Make** utility installed [directly](https://gnuwin32.sourceforge.net/packages/make.htm) or via [conda](https://anaconda.org/conda-forge/make)

### Optional

- **NVIDIA GPU Computing Toolkit** [directly](https://developer.nvidia.com/cuda-downloads) or via [conda](https://anaconda.org/nvidia) | Expected location:

  [ Windows ] `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\YOUR_VERSION\lib\x64\`

  [ Linux ]  `$CONDA_PREFIX/lib`

- [CVODE](https://computing.llnl.gov/projects/sundials/sundials-software) | Expected location:

  [ Windows ] `C:\Program Files (x86)\sundials-4.1.0\instdir`

  [ Linux ] `/usr/local/sundials-4.1.0/instdir`

  A guide for custom installation of sundials-4.1.0 from this repository can be found [here](#install-cvode-from-sundials-4.1.0) 

  **Note:** Shared libraries for Fortran modules are not distributed in conda package of [sundials](https://anaconda.org/conda-forge/sundials). CVODE has to be built beforehand locally when linked dynamically. During runtime, path has to added to LD_LIBRARY_PATH on Linux.

### Compilation

Prepare your terminal, so that ifort compiler can be found:

- [ Windows]

  ```bash
  "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
  ```

- [ Linux ]

  ```bash
  . /opt/intel/oneapi/setvars.sh
  ```

Set flags and paths in [Makefile](https://github.com/cmt-dtu-energy/MagTense/blob/master/Makefile) corresponding to your setup and then run:

```bash
cd path/to/MagTense/
make USE_CUDA=1 USE_CVODE=1 USE_MATLAB=1 CVODE_ROOT=/usr/local/sundials-4.1.0/instdir MATLAB_INCLUDE=/usr/local/MATLAB/R2021b/extern/include MKL_ROOT=/opt/intel/oneapi/mkl/latest     
```

## Building MEX-files

Prepare your terminal, so that ifort compiler can be found:

- [ Windows]

  ```bash
  "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
  ```

  Modify variable `CMDLINE100` in `C:\Users\<user>\AppData\Roaming\MathWorks\MATLAB\R2023b\mex_FORTRAN_win64.xml`, so that the paths stored in `INCLUDE` can be found during building.
  
  ```bash
  CMDLINE100="$COMPILER /c $COMPFLAGS $INCLUDE $OPTIM $SRC /Fo$OBJ"
  ```

  Configuration of FORTRAN compiler used for mex:

  ```bash
  mex -setup fortran
  ```

- [ Linux ]

  ```bash
  . /opt/intel/oneapi/setvars.sh
  ```

Start Matlab via command line and run build script:

```bash
cd path/to/MagTense/matlab
/path/to/Matlab/installation/bin/matlab -nodisplay -nosplash -nodesktop
run('buildMagTenseMEX.m')
```

## Runtime

Prepare your terminal, so that ifort compiler and required libraries can be found:

- [ Windows]

  ```bash
  "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
  ```

- [ Linux ]

  ```bash
  . /opt/intel/oneapi/setvars.sh
  export LD_LIBRARY_PATH=<path_to_CVODE_libs>:<path_to_CUDA_libs>:$LD_LIBRARY_PATH
  ```

- Test with standard problem #3

  ```bash
  cd path/to/MagTense/matlab/
  /path/to/Matlab/installation/bin/matlab -nodisplay -nosplash -nodesktop
  run('examples/Micromagnetism/mumag_micromag_Std_problem_3/Standard_problem_3.m')
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

## Compile changes in CUDA code on Windows with Intel Fortran

Intel's Fortran Compiler and the PG Fortran compiler (which is used for doing CUDA directly in Fortran) are not compatible. One simply cannot link object files from the two compilers as there is no standard for this.

The strategy then is to make the CUDA GPU kernel in C++, compile this with the nvcc compiler (that uses MSVC at the core), output to an object file and then compile a C++ wrapper with Intel C++ compiler (icx) that includes and uses the nvcc compiled file. The output here should be an `obj` file that can be called from Fortran (using the Intel Fortran compiler) via the standard way (iso_c_binding) for calling C++ functions from Fortran.

The following details the required steps if changes has been made to the CUDA code.

### Step 0

- Install VS 2019 with Intel's C++ and Fortran compilers (and MKL as this is used also by MagTense) and with MSVC compiler
- Install [CUDA SDK](https://developer.nvidia.com/cuda-downloads)
- Check that the PATH variable is pointing to the latest version of CUDA. This can be set using the registry at Computer\HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment\Path

### Step 1

From a standard command prompt, navigate to the MagTense dir "MagTense/source/MagTenseFortranCUDA/cuda" and compile with nvcc, using a specific CUDA version and the latest version of Visual Studio: 

```bash
"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin\nvcc.exe" -c MagTenseCudaBlas.cu -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.36.32532\bin\Hostx64\x64"
```

### Step 2

From the Intel 64 prompt (Start Menu -> Intel oneAPI 20XX -> Intel oneAPI command prompt for Intel 64 for Visual Studio XX) compile the C++ wrapper with icx including the cuda stuff:

```bash
icx -c MagTenseCudaBlasICLWrapper.cxx
```

### Step 3

Compile MagTense as usual. There is a project in the MagTense solution that deals with CUDA.

### Step 4. (optional)

Check the dependencies of the compiled MEX-file using the Visual Studio dumpbin utility:

```bash
"C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.36.32532\bin\Hostx64\x64\dumpbin.exe" /dependents MagTenseLandauLifshitzSolver_mex.mexw64
```
