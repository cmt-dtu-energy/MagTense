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
ForceIntegrator_path = '../source/MagneticForceIntegrator/MagneticForceIntegrator';
FortranCuda_path = '../source/MagTenseFortranCuda/cuda';

if (ispc)
    VS_STUDIO = true;
    MKL_STATIC = true;
    mkl_include = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include"';
    mkl_lp64 = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\intel64\lp64"';
    mkl_lib = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib\intel64"';
    cuda_root = '"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\lib\x64"';
    cvode_root = '"C:\Program Files (x86)\sundials-4.1.0\instdir"';
    mex_suffix = 'w';
else
    compiler_root = '/opt/intel/oneapi/compiler/latest/linux';
    mkl_root = '/opt/intel/oneapi/mkl/latest';
    cuda_root = join([getenv('CONDA_PREFIX') '/lib'], '');
    cvode_root = '/usr/local/sundials-4.1.0/instdir';
    mex_suffix = 'a';
end

if (USE_RELEASE)
    DEBUG = '-v';
    BUILD = '/x64/Release';
else
    DEBUG = '-v -g';
    BUILD = '/x64/Debug';
end

BUILD_NO_CUDA_MagTenseMicroMag = [BUILD '_no_CUDA'];
if (USE_CUDA)
    CUDA = ['-L' cuda_root ' -lcublas -lcudart -lcuda -lcusparse'];
    if (ispc)
        OBJS = ['OBJS="$OBJS ' FortranCuda_path '/MagTenseCudaBlasICLWrapper.obj ' FortranCuda_path '/MagTenseCudaBlas.obj" '];
    else
        OBJS = ['OBJS="$OBJS ' FortranCuda_path '/MagTenseCudaBlasICLWrapper.o' FortranCuda_path '/MagTenseCudaBlas.o" '];
    end
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

if (VS_STUDIO)
    BUILD = [BUILD ' '];
    BUILD_MagTenseMicroMag = [BUILD_MagTenseMicroMag ' '];
    BUILD_NO_CUDA_MagTenseMicroMag = [BUILD_NO_CUDA_MagTenseMicroMag ' '];
else
    BUILD = ' ';
    BUILD_MagTenseMicroMag = ' ';
    BUILD_NO_CUDA_MagTenseMicroMag = ' ';
end


if (ispc)
    DEFINES = '-R2018a';
    FFLAGS = 'COMPFLAGS="$COMPFLAGS /free /nologo /real-size:64 /O2 /assume:nocc_omp /Qopenmp /fpp /fpe:0 /fp:source /fp:precise /libs:static"';
    INCLUDE = ['-I' mkl_include ' -I' mkl_lp64 ' -I' NumericalIntegration_path BUILD '-I' DemagField_path BUILD '-I' TileDemagTensor_path ...
        BUILD '-I' MagTenseMicroMag_path BUILD_MagTenseMicroMag '-I' ForceIntegrator_path BUILD CVODE_include];
    LIBS = ['-L' MagTenseMicroMag_path BUILD_MagTenseMicroMag '-lMagTenseMicroMag -L' DemagField_path BUILD ...
        ' -lDemagField -L' TileDemagTensor_path BUILD ' -lTileDemagTensor -L' NumericalIntegration_path BUILD ...
        ' -lNumericalIntegration -L' ForceIntegrator_path BUILD ' -lMagneticForceIntegrator'];
    INCLUDE_NO_CUDA = ['-I' mkl_include ' -I' mkl_lp64 ' -I' NumericalIntegration_path BUILD '-I' DemagField_path BUILD '-I' TileDemagTensor_path ...
        BUILD '-I' MagTenseMicroMag_path BUILD_NO_CUDA_MagTenseMicroMag '-I' ForceIntegrator_path BUILD CVODE_include];
    LIBS_NO_CUDA = ['-L' MagTenseMicroMag_path BUILD_NO_CUDA_MagTenseMicroMag ' -lMagTenseMicroMag -L' DemagField_path BUILD ...
        ' -lDemagField -L' TileDemagTensor_path BUILD ' -lTileDemagTensor -L' NumericalIntegration_path BUILD ...
        ' -lNumericalIntegration -L' ForceIntegrator_path BUILD ' -lMagneticForceIntegrator'];
    
    if (MKL_STATIC)
        MKL = ['-L' mkl_lib ' -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -lmkl_blas95_lp64 -llibiomp5md'];
    else
        MKL = ['-L' mkl_lib ' -lmkl_rt -lmkl_blas95_lp64'];
    end
else
    DEFINES = ['FC="' compiler_root '/bin/intel64/ifort" DEFINES="-DMATLAB_DEFAULT_RELEASE=R2018a"'];
    FFLAGS = 'FFLAGS="-r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source -fp-model precise -fpic -libs:static"';
    INCLUDE = ['INCLUDE="$INCLUDE -I' mkl_root '/include -I' mkl_root '/include/intel64/lp64 -I' ... 
        NumericalIntegration_path ' -I' DemagField_path ' -I' TileDemagTensor_path ' -I' MagTenseMicroMag_path ...
        ' -I' ForceIntegrator_path ' ' CVODE_include '"'];
    LIBS = ['-L' MagTenseMicroMag_path ' -lMagTenseMicroMag -L' DemagField_path ' -lDemagField -L' ...
        TileDemagTensor_path ' -lTileDemagTensor -L' NumericalIntegration_path ' -lNumericalIntegration -L' ...
        ForceIntegrator_path ' -lMagneticForceIntegrator'];
    MKL = ['-L' mkl_root '/lib/intel64 -lmkl_rt -lmkl_blas95_lp64 -lpthread -lm -ldl'];
end

%%------------------------------------------------------------------
%%--------------- Build the MEX files ------------------------------
%%----------------------------------- ------------------------------
if (ispc) && (VS_STUDIO)
    names = ["MagTenseLandauLifshitzSolver", "MagTenseLandauLifshitzSolverNoCUDA", "IterateMagnetization", "getHFromTiles", "getNFromTile", "getMagForce"];
else
    names = ["MagTenseLandauLifshitzSolver", "IterateMagnetization", "getHFromTiles", "getNFromTile", "getMagForce"];
end

for i = 1:length(names)
    if names(i) == "MagTenseLandauLifshitzSolverNoCUDA"
        source = [mex_root 'MagTenseLandauLifshitzSolver_mex.f90'];
        mex_str = ['mex' ' ' DEBUG ' ' DEFINES ' ' FFLAGS ' ' INCLUDE_NO_CUDA ' ' OBJS ' ' LIBS_NO_CUDA ' ' MKL ' ' CUDA ' ' CVODE ' ' join(source, '')];
        orig_name = "MagTenseLandauLifshitzSolver";
    else
        source = [mex_root names(i) '_mex.f90'];
        mex_str = ['mex' DEBUG DEFINES FFLAGS INCLUDE OBJS LIBS MKL CUDA CVODE join(source, '')];
        orig_name = names(i);
    end
    disp(join(mex_str, ' '))
    eval_MEX(join(mex_str, ' '))
    pause(pause_time)
    movefile(join([orig_name '_mex.mex' mex_suffix '64'], ''), join(['MEX_files/' names(i) '_mex.mex' mex_suffix '64'], ''));
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
