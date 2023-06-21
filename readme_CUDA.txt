How to compile changes in the CUDA code in Intel Fortran in MagTense on Windows:

Intel's Fortran Compiler and the PG Fortran compiler (which is used for doing CUDA directly in Fortran) are not compatible. One simply cannot link object files from the two compilers as there is no standard for this.

The strategy then is to make the CUDA GPU kernel in C++, compile this with the nvcc compiler (that uses MSVC at the core), output to an object file and then compile a C++ wrapper with Intel C++ compiler (icx) that includes and uses the nvcc compiled file. The output here should be an obj file that can be called from Fortran (using the Intel Fortran compiler) via the standard way (iso_c_binding) for calling C++ functions from Fortran.

The following details the required steps if changes has been made to the CUDA code.
 
Step 0.
 Install VS 2019 with Intel's C++ and Fortran compilers (and MKL as this is used also by MagTense) and with MSVC compiler
 Install CUDA SDK: https://developer.nvidia.com/cuda-downloads
 Check that the PATH variable is pointing to the latest version of CUDA. This can be set using the registry at Computer\HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment\Path

Step 1. 
 From a standard command prompt, navigate to the MagTense dir "MagTense/source/MagTenseFortranCUDA/cuda" and compile with nvcc, using a specific CUDA version and the latest version of Visual Studio: 
 
 "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin\nvcc.exe" -c MagTenseCudaBlas.cu -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.36.32532\bin\Hostx64\x64"
 
Step 2.
 From the Intel 64 prompt (Start Menu -> Intel oneAPI 20XX -> Intel oneAPI command prombt for Intel 64 for Visual Studio XX compile the C++ wrapper with icx including the cuda stuff*:
 
 icx -c MagTenseCudaBlasICLWrapper.cxx
 
Step 3.
 Compile MagTense as usual (there is a project in the MagTense solution that deals with CUDA)
 
Step 4. (optional)
 Check the dependencies of the compiled MEX-file using the Visual Studio dumpbin utility:

 "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.36.32532\bin\Hostx64\x64\dumpbin.exe" /dependents MagTenseLandauLifshitzSolver_mex.mexw64