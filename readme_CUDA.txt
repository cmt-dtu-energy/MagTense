How to compile changes in the CUDA code in Intel Fortran in MagTense on Windows:

Intel's Fortran Compiler and the PG Fortran compiler (which is used for doing CUDA directly in Fortran) are not compatible. One simply cannot link object files from the two compilers as there is no standard for this.

The strategy then is to make the CUDA GPU kernel in C++, compile this with the nvcc compiler (that uses MSVC at the core), output to an object file and then compile a C++ wrapper with Intel C++ compiler (icl) that includes and uses the nvcc compiled file. The output here should be an obj file that can be called from Fortran (using the Intel Fortran compiler) via the standard way (iso_c_binding) for calling C++ functions from Fortran.

The following details the required steps if changes has been made to the CUDA code.
 
Step 0.
 Install VS 2019 with Intel's C++ and Fortran compilers (and MKL as this is used also by MagTense) and with MSVC compiler
 Install CUDA SDK: https://developer.nvidia.com/cuda-downloads

Step 1. 
 From the MS Visual Studio x64 native prompt navigate to the MagTense dir and compile with nvcc: 
 
 cd MagTense/source/MagTenseFortranCUDA/cuda
 nvcc -c MagTenseCudaBlas.cu
 
Step 2.
 From the Intel 64 prompt (Start Menu -> Intel Parallel Studio -> latest x64 cmd) compile the C++ wrapper with icl including the cuda stuff:
 
 icl -c MagTenseCudaBlasICLWrapper.cxx MagTenseCudaBlas.obj  /link "c:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.2\lib\x64\cublas.lib"
 
Step 3.
 Compile MagTense as usual (there is a project in the MagTense solution that deals with CUDA)
 
 