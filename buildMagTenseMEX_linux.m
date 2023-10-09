function buildMagTenseMEX()
%use clear all as this also clears dependencies to the .mex files and thus they can be overwritten


%{
## mex_FORTRAN_glnxa64"
/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_blas95_lp64.a -Wl,--start-group 
/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.a 
/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.a 
/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.a -Wl,--end-group

### NO CUDA
### DYNAMIC MATLAB -- WORKS
mex -v FC="ifort" 
DEFINES="-DMATLAB_DEFAULT_RELEASE=R2018a -DUSE_MEX_CMD -DMKL_LP64" 
FFLAGS="-r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source -fp-model precise -fpic -libs:static" 
INCLUDE="$INCLUDE -I/opt/intel/oneapi/mkl/latest/include -I/opt/intel/oneapi/mkl/latest/include/intel64/ilp64 
-Isource/NumericalIntegration/NumericalIntegration -Isource/TileDemagTensor/TileDemagTensor -Isource/DemagField/DemagField -Isource/MagTenseMicroMag" 
-Lsource/MagTenseMicroMag/ -lMagTenseMicroMag -Lsource/DemagField/DemagField -lDemagField -Lsource/TileDemagTensor/TileDemagTensor -lTileDemagTensor 
-Lsource/NumericalIntegration/NumericalIntegration -lNumericalIntegration -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64 -lpthread -lm -ldl -Lsource/MagTenseMicroMag/ -lMagTenseMicroMag 
source/MagTenseMEX/MagTenseMEX/MagTenseLandauLifshitzSolver_mex.f90


### WITH CUDA
mex -v FC="ifort" 
DEFINES="-DMATLAB_DEFAULT_RELEASE=R2018a -DUSE_MEX_CMD -DMKL_LP64" 
FFLAGS="-r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source -fp-model precise -fpic -libs:static" 
INCLUDE="$INCLUDE -I/opt/intel/oneapi/mkl/latest/include -I/opt/intel/oneapi/mkl/latest/include/intel64/ilp64 
-Isource/NumericalIntegration/NumericalIntegration -Isource/TileDemagTensor/TileDemagTensor -Isource/DemagField/DemagField -Isource/MagTenseMicroMag"
OBJS="$OBJS source/MagTenseFortranCuda/cuda/MagTenseCudaBlasICLWrapper.o source/MagTenseFortranCuda/cuda/MagTenseCudaBlas.o"
LINKLIBS="$LINKLIBS -L/opt/intel/oneapi/mkl/latest/lib/intel64 -L/opt/intel/oneapi/compiler/latest/lib/intel64_lin -L' cuda_root '/home/spol/miniconda3/envs/mt-cuda-py310/lib"'
-Lsource/MagTenseMicroMag/ -lMagTenseMicroMag -Lsource/DemagField/DemagField -lDemagField -Lsource/TileDemagTensor/TileDemagTensor -lTileDemagTensor 
-Lsource/NumericalIntegration/NumericalIntegration -lNumericalIntegration -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64 -lpthread -lm -ldl 
-L/home/spol/miniconda3/envs/mt-py311-cuda12/lib -lcublas -lcudart -lcuda -lcusparse -Lsource/MagTenseMicroMag/ -lMagTenseMicroMag 
source/MagTenseMEX/MagTenseMEX/MagTenseLandauLifshitzSolver_mex.f90 
%}

compiler_root = '/opt/intel/oneapi/compiler/latest';
mkl_root = '/opt/intel/oneapi/mkl/latest';
cuda_root = '/home/spol/miniconda3/envs/mt-cuda-py310';
% cvode_root = '/home/spol/cvode-6.6.0/instdir';

mex_root = 'source/MagTenseMEX/MagTenseMEX/';
USE_CUDA = 0;
RELEASE = 1;
% USE_CVODE = 0;
mex_suffix = 'a';
pause_time = 1;

if (RELEASE)
    DEBUG = '-v';
    build_str = 'linux/Release';
else
    DEBUG = '-v -g';
    build_str = 'linux/Debug';
end

if (USE_CUDA)
    CUDA_FLAG = '-DUSE_CUDA=1';
    CUDA = ['-L' cuda_root '/lib -lcublas -lcudart -lcuda -lcusparse'];
    OBJS = 'OBJS="$OBJS source/MagTenseFortranCuda/cuda/MagTenseCudaBlasICLWrapper.o source/MagTenseFortranCuda/cuda/MagTenseCudaBlas.o" ';
    name_suffix = '_mex';
else
    CUDA_FLAG = '-DUSE_CUDA=0';
    CUDA = '';
    OBJS = '';
    name_suffix = 'NoCUDA_mex';
end

NumericalIntegration_include = '-Isource/NumericalIntegration/NumericalIntegration';
NumericalIntegration_lib = '-Lsource/NumericalIntegration/NumericalIntegration -lNumericalIntegration';

DemagField_include = '-Isource/DemagField/DemagField';
DemagField_lib = '-Lsource/DemagField/DemagField -lDemagField';

TileDemagTensor_include = '-Isource/TileDemagTensor/TileDemagTensor';
TileDemagTensor_lib = '-Lsource/TileDemagTensor/TileDemagTensor -lTileDemagTensor';

MagTenseMicroMag_include = '-Isource/MagTenseMicroMag';
MagTenseMicroMag_lib = '-Lsource/MagTenseMicroMag/ -lMagTenseMicroMag';

COMPILE = ['FC="' compiler_root '/linux/bin/intel64/ifort"'];

DEFINES = 'DEFINES="-DMATLAB_DEFAULT_RELEASE=R2018a -DUSE_MEX_CMD -DMKL_LP64"';

FFLAGS = 'FFLAGS="-i8 -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source -fp-model precise -fpic -libs:static"';

INCLUDE = ['INCLUDE="$INCLUDE -I' mkl_root '/include -I' mkl_root '/include/intel64/ilp64 ' ... 
NumericalIntegration_include ' ' DemagField_include ' ' TileDemagTensor_include ' ' MagTenseMicroMag_include '"'];

LINKLIBS = ['LINKLIBS="$LINKLIBS -L' mkl_root '/lib/intel64 -L' compiler_root '/lib/intel64_lin -L' cuda_root '/lib"'];

MKL = ['-L' mkl_root '/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64 -lpthread -lm -ldl'];

%{
 if (USE_CVODE)
    CVODE_str = [cvode_root '/lib -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_'...
    'sunlinsoldense -lsundials_fnvecserial_mod -lsundials_cvode -lsundials_fsunnonlinsolfixedpoint_mod -I' cvode_root '/include'];
else
    CVODE_str = '';
    build_str = [build_str  '_no_CVODE'];
    build_str_MagTenseMicroMag = [build_str_MagTenseMicroMag '_no_CVODE'];
end 
%}

% MagneticForceIntegrator_str = ['-Lsource/MagneticForceIntegrator/MagneticForceIntegrator/' build_str '/ -lMagneticForceIntegrator '...
% '-Lsource/MagneticForceIntegrator/MagneticForceIntegrator/MagneticForceIntegratorLib.a -Isource/MagneticForceIntegrator/MagneticForceIntegrator/' build_str '/'];


%%------------------------------------------------------------------
%%--------------- Build the MEX files ------------------------------
%%----------------------------------- ------------------------------
name = 'MagTenseLandauLifshitzSolver';
source = [mex_root name '_mex.f90'];

mex_str = ['mex ' DEBUG ' ' COMPILE ' ' DEFINES ' ' FFLAGS ' ' INCLUDE ' ' OBJS LINKLIBS ' ' ...
MagTenseMicroMag_lib ' ' DemagField_lib ' ' TileDemagTensor_lib ' ' NumericalIntegration_lib ' ' MKL ' ' CUDA ' ' source];

eval(join(mex_str, ' '))
pause(pause_time)
movefile([name '_mex.mex' mex_suffix '64'], ['matlab/MEX_files/' name name_suffix '.mex' mex_suffix '64']);

%{
source = 'IterateMagnetization_mex';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' DemagField_str ' ' TileDemagTensor_str ' ' NumericalIntegration_str ' ' mex_root source '.f90 ' Options_str];
eval(mex_str)
pause(pause_time)
movefile(['IterateMagnetization_mex.mex' MEX_str '64'],['matlab/MEX_files/IterateMagnetization_mex.mex' MEX_str '64']);
        
source = 'getHFromTiles_mex';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' DemagField_str ' ' TileDemagTensor_str ' ' NumericalIntegration_str ' ' mex_root source '.f90 ' Options_str];
eval(mex_str)
pause(pause_time)
movefile(['getHFromTiles_mex.mex' MEX_str '64'],['matlab/MEX_files/getHFromTiles_mex.mex' MEX_str '64']);
    
source = 'getNFromTile_mex';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' mex_root source '.f90 ' Options_str];
eval(mex_str)
pause(pause_time)
movefile(['getNFromTile_mex.mex' MEX_str '64'],['matlab/MEX_files/getNFromTile_mex.mex' MEX_str '64']);
    
source = 'getMagForce_mex';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' MagneticForceIntegrator_str ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' mex_root source '.f90 ' Options_str];
eval(mex_str)
pause(pause_time)
movefile(['getMagForce_mex.mex' MEX_str '64'],['matlab/MEX_files/getMagForce_mex.mex' MEX_str '64']);
    
%}
 
end
