function buildMagTenseMEX(USE_RELEASE, USE_CUDA, USE_CVODE)
%use clear all as this also clears dependencies to the .mex files and thus they can be overwritten

arguments
    USE_RELEASE {mustBeNumericOrLogical} = true;
    USE_CUDA {mustBeNumericOrLogical} = true;
    USE_CVODE {mustBeNumericOrLogical} = true;
end

compiler_root = '/opt/intel/oneapi/compiler/latest';
mkl_root = '/opt/intel/oneapi/mkl/latest';
cuda_root = '/home/spol/miniconda3/envs/mt-cuda-py310';
cvode_root = '/home/spol/cvode-6.6.0/instdir';
mex_root = 'source/MagTenseMEX/MagTenseMEX/';

NumericalIntegration_path = 'source/NumericalIntegration/NumericalIntegration';
DemagField_path = 'source/DemagField/DemagField';
TileDemagTensor_path = 'source/TileDemagTensor/TileDemagTensor';
MagTenseMicroMag_path = 'source/MagTenseMicroMag';
ForceIntegrator_path = 'source/MagneticForceIntegrator/MagneticForceIntegrator/';

mex_suffix = 'a';

if (USE_RELEASE)
    DEBUG = '-v';
else
    DEBUG = '-v -g';
end

name_suffix = '';
if (USE_CUDA)
    CUDA = ['-L' cuda_root '/lib -lcublas -lcudart -lcuda -lcusparse'];
    OBJS = 'OBJS="$OBJS source/MagTenseFortranCuda/cuda/MagTenseCudaBlasICLWrapper.o source/MagTenseFortranCuda/cuda/MagTenseCudaBlas.o" ';
else
    CUDA = '';
    OBJS = '';
    name_suffix = 'NoCUDA';
end

%{
 if (USE_CVODE)
    CVODE_include = ['-I' cvode_root '/include'];
    CVODE = ['-L' cvode_root '/lib -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense' ...
    ' -lsundials_fnvecserial_mod -lsundials_cvode -lsundials_fsunnonlinsolfixedpoint_mod'];
else
    CVODE = '';
    name_suffix = [name_suffix '_no_CVODE'];
end 
%}

CVODE = '';
COMPILE = ['FC="' compiler_root '/linux/bin/intel64/ifort"'];
DEFINES = 'DEFINES="-DMATLAB_DEFAULT_RELEASE=R2018a"';
FFLAGS = 'FFLAGS="-i8 -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source -fp-model precise -fpic -libs:static"';
INCLUDE = ['INCLUDE="$INCLUDE -I' mkl_root '/include -I' mkl_root '/include/intel64/ilp64 -I' ... 
    NumericalIntegration_path ' -I' DemagField_path ' -I' TileDemagTensor_path ' -I' MagTenseMicroMag_path ...
    ' -I' ForceIntegrator_path '"'];
LIBS = ['-L' MagTenseMicroMag_path ' -lMagTenseMicroMag -L' DemagField_path ' -lDemagField -L' ...
    TileDemagTensor_path ' -lTileDemagTensor -L' NumericalIntegration_path ' -lNumericalIntegration -L' ...
    ForceIntegrator_path ' -lMagneticForceIntegrator'];
MKL = ['-L' mkl_root '/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64 -lpthread -lm -ldl'];

%%------------------------------------------------------------------
%%--------------- Build the MEX files ------------------------------
%%----------------------------------- ------------------------------
names = ["MagTenseLandauLifshitzSolver", "IterateMagnetization", "getHFromTiles", "getNFromTile", "getMagForce"];

for i = 1:length(names)
    source = [mex_root names(i) '_mex.f90'];
    mex_str = ['mex' DEBUG COMPILE DEFINES FFLAGS INCLUDE OBJS LIBS MKL CUDA CVODE join(source, '')];
    disp(join(mex_str, ' '))
    eval(join(mex_str, ' '))
    if names(i) ~= "MagTenseLandauLifshitzSolver"
        name_suffix = '';
    end
    movefile(join([names(i) '_mex.mex' mex_suffix '64'], ''), join(['matlab/MEX_files/' names(i) name_suffix '_mex.mex' mex_suffix '64'], ''));
end
 
end
