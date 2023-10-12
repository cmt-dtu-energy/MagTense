function buildMagTenseMEX(USE_RELEASE, USE_CUDA, USE_CVODE)
%use clear all as this also clears dependencies to the .mex files and thus they can be overwritten

arguments
    USE_RELEASE {mustBeNumericOrLogical} = true;
    USE_CUDA {mustBeNumericOrLogical} = true;
    USE_CVODE {mustBeNumericOrLogical} = true;
end

pause_time = 1; %Time to wait between making and moving the generated files
mex_root = '../source/MagTenseMEX/MagTenseMEX/';
NumericalIntegration_path = '../source/NumericalIntegration/NumericalIntegration';
DemagField_path = '../source/DemagField/DemagField';
TileDemagTensor_path = '../source/TileDemagTensor/TileDemagTensor';
MagTenseMicroMag_path = '../source/MagTenseMicroMag';
ForceIntegrator_path = '../source/MagneticForceIntegrator/MagneticForceIntegrator/';

if (ispc)
    compiler_root = '"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows"';
    mkl_root = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest"';
    cuda_root = '"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\lib\x64"';
    cvode_root = '"C:\Program Files (x86)\sundials-4.1.0\instdir"';
    mex_suffix = 'w';
    setenv('IFORT_COMPILER19', compiler_root); % Necessary for Matlab to recognize oneAPI
else
    compiler_root = '/opt/intel/oneapi/compiler/latest/linux';
    mkl_root = '/opt/intel/oneapi/mkl/latest';
    cuda_root = join([getenv('CONDA_PREFIX') '/lib'], '');
    cvode_root = '/usr/local/sundials-4.1.0/instdir';
    mex_suffix = 'a';
end

if (USE_RELEASE)
    DEBUG = '-v';
    BUILD = 'x64/Release';
else
    DEBUG = '-v -g';
    BUILD = 'x64/Debug';
end

BUILD_NO_CUDA_MagTenseMicroMag = [BUILD '_no_CUDA'];
if (USE_CUDA)
    CUDA = ['-L' cuda_root ' -lcublas -lcudart -lcuda -lcusparse'];
    OBJS = 'OBJS="$OBJS ../source/MagTenseFortranCuda/cuda/MagTenseCudaBlasICLWrapper.o ../source/MagTenseFortranCuda/cuda/MagTenseCudaBlas.o" ';
    BUILD_MagTenseMicroMag = BUILD;
else
    CUDA = '';
    OBJS = '';
    BUILD_MagTenseMicroMag = [BUILD '_no_CUDA'];
end

if (USE_CVODE)
    CVODE_include = join(['-I' cvode_root '/fortran'], '');
    CVODE = ['-L' cvode_root '/lib -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense' ...
    ' -lsundials_fnvecserial_mod -lsundials_cvode -lsundials_fsunnonlinsolfixedpoint_mod'];
else
    CVODE_include = '';
    CVODE = '';
    BUILD = [BUILD '_no_CVODE'];
    BUILD_MagTenseMicroMag = [BUILD_MagTenseMicroMag '_no_CVODE'];
    BUILD_NO_CUDA_MagTenseMicroMag = [BUILD_NO_CUDA_MagTenseMicroMag '_no_CVODE'];
end

COMPILE = ['FC="' compiler_root '/bin/intel64/ifort"'];
DEFINES = 'DEFINES="-DMATLAB_DEFAULT_RELEASE=R2018a"';
FFLAGS = 'FFLAGS="-i8 -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source -fp-model precise -fpic -libs:static"';

if (ispc)
    INCLUDE = ['INCLUDE="$INCLUDE -I' mkl_root '/include -I' mkl_root '/include/intel64/ilp64 -I' ... 
        NumericalIntegration_path '/' BUILD ' -I' DemagField_path '/' BUILD ' -I' TileDemagTensor_path '/' BUILD ...
        ' -I' MagTenseMicroMag_path '/' BUILD_MagTenseMicroMag ' -I' ForceIntegrator_path '/' BUILD ' ' CVODE_include '"'];
    LIBS = ['-L' MagTenseMicroMag_path '/' BUILD_MagTenseMicroMag ' -lMagTenseMicroMag -L' DemagField_path '/' BUILD ...
        ' -lDemagField -L' TileDemagTensor_path '/' BUILD ' -lTileDemagTensor -L' NumericalIntegration_path '/' BUILD ...
        ' -lNumericalIntegration -L' ForceIntegrator_path '/' BUILD ' -lMagneticForceIntegrator'];
    INCLUDE_NO_CUDA = ['INCLUDE="$INCLUDE -I' mkl_root '/include -I' mkl_root '/include/intel64/ilp64 -I' ... 
        NumericalIntegration_path '/' BUILD ' -I' DemagField_path '/' BUILD ' -I' TileDemagTensor_path '/' BUILD ...
        ' -I' MagTenseMicroMag_path '/' BUILD_NO_CUDA_MagTenseMicroMag ' -I' ForceIntegrator_path '/' BUILD ' ' CVODE_include '"'];
    LIBS_NO_CUDA = ['-L' MagTenseMicroMag_path '/' BUILD_NO_CUDA_MagTenseMicroMag ' -lMagTenseMicroMag -L' DemagField_path '/' BUILD ...
        ' -lDemagField -L' TileDemagTensor_path '/' BUILD ' -lTileDemagTensor -L' NumericalIntegration_path '/' BUILD ...
        ' -lNumericalIntegration -L' ForceIntegrator_path '/' BUILD ' -lMagneticForceIntegrator'];
else
    INCLUDE = ['INCLUDE="$INCLUDE -I' mkl_root '/include -I' mkl_root '/include/intel64/ilp64 -I' ... 
        NumericalIntegration_path ' -I' DemagField_path ' -I' TileDemagTensor_path ' -I' MagTenseMicroMag_path ...
        ' -I' ForceIntegrator_path ' ' CVODE_include '"'];
    LIBS = ['-L' MagTenseMicroMag_path ' -lMagTenseMicroMag -L' DemagField_path ' -lDemagField -L' ...
        TileDemagTensor_path ' -lTileDemagTensor -L' NumericalIntegration_path ' -lNumericalIntegration -L' ...
        ForceIntegrator_path ' -lMagneticForceIntegrator'];
end

MKL = ['-L' mkl_root '/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64 -lpthread -lm -ldl'];

%%------------------------------------------------------------------
%%--------------- Build the MEX files ------------------------------
%%----------------------------------- ------------------------------
if (ispc)
    names = ["MagTenseLandauLifshitzSolver", "MagTenseLandauLifshitzSolverNoCUDA", "IterateMagnetization", "getHFromTiles", "getNFromTile", "getMagForce"];
else
    names = ["MagTenseLandauLifshitzSolver", "IterateMagnetization", "getHFromTiles", "getNFromTile", "getMagForce"];
end

for i = 1:length(names)
    if names(i) == "MagTenseLandauLifshitzSolverNoCUDA"
        source = [mex_root 'MagTenseLandauLifshitzSolver_mex.f90'];
        mex_str = ['mex' DEBUG COMPILE DEFINES FFLAGS INCLUDE_NO_CUDA OBJS LIBS_NO_CUDA MKL CUDA CVODE join(source, '')];
    else
        source = [mex_root names(i) '_mex.f90'];
        mex_str = ['mex' DEBUG COMPILE DEFINES FFLAGS INCLUDE OBJS LIBS MKL CUDA CVODE join(source, '')];
    end
    disp(join(mex_str, ' '))
    eval_MEX(join(mex_str, ' '))
    pause(pause_time)
    movefile(join([names(i) '_mex.mex' mex_suffix '64'], ''), join(['matlab/MEX_files/' names(i) '_mex.mex' mex_suffix '64'], ''));
end
 
end

function eval_MEX(mex_str)
try
    eval(mex_str);
catch ME
    if (strcmp(ME.message(91:117),'mt : general error c101008d'))
        fail_mex = true; 
        while fail_mex
            try 
                disp('Microsoft manifest tool error - retrying')
                eval(mex_str); 
                fail_mex = false; 
            catch
                continue
            end
        end
    else
        disp(ME.message)
        rethrow(ME)
    end
end
end
