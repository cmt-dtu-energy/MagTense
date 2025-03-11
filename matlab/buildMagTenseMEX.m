function buildMagTenseMEX(USE_RELEASE, USE_CUDA, USE_CVODE)
%use clear all as this also clears dependencies to the .mex files and thus they can be overwritten

arguments
    USE_RELEASE {mustBeNumericOrLogical} = true;
    USE_CUDA {mustBeNumericOrLogical} = false;
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
    mkl_lp64 = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\mkl\intel64\lp64"';
    mkl_lib = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib"';
    cuda_root = '"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\lib\x64"';
    cvode_include = '"C:\Program Files (x86)\sundials-4.1.0\instdir\fortran"';
    cvode_lib = '"C:\Program Files (x86)\sundials-4.1.0\instdir\lib"';
    mex_suffix = 'w';
else
    VS_STUDIO = false;
    MKL_STATIC = true;
    compiler_root = '/usr/share/miniconda/envs/magtense-env';
    mkl_root = '/usr/share/miniconda/envs/magtense-env';
    mkl_lib = '/usr/share/miniconda/envs/magtense-env/lib';
    cuda_root = '/home/spol/miniconda3/envs/magtense-env-py12/lib';
    cvode_include = '/usr/local/sundials-4.1.0/instdir/fortran';
    cvode_lib = '/usr/local/sundials-4.1.0/instdir/lib';
    mex_suffix = 'a';
end

if (USE_RELEASE)
    DEBUG = '';
    BUILD = '/x64/Release';
else
    DEBUG = '-g';
    BUILD = '/x64/Debug';
end

if (USE_CUDA)
    CUDA = ['-L' cuda_root ' -lcublas -lcudart -lcuda -lcusparse'];
    if (ispc)
        OBJS = ['OBJS="$OBJS ' FortranCuda_path '/MagTenseCudaBlasICLWrapper.obj ' FortranCuda_path '/MagTenseCudaBlas.obj" '];
    else
        OBJS = ['OBJS="$OBJS ' FortranCuda_path '/MagTenseCudaBlasICLWrapper.o ' FortranCuda_path '/MagTenseCudaBlas.o" '];
    end
    BUILD_MagTenseMicroMag = BUILD;
else
    CUDA = '';
    OBJS = '';
    BUILD_MagTenseMicroMag = [BUILD '_no_CUDA'];
end

if (USE_CVODE)
    CVODE_include = join(['-I' cvode_include], '');
    CVODE = ['LINKLIBS="$LINKLIBS ' cvode_lib '/libsundials_nvecserial.a ' cvode_lib '/libsundials_sunmatrixdense.a ' ...
        cvode_lib '/libsundials_sunlinsoldense.a ' cvode_lib '/libsundials_fnvecserial_mod.a ' ...
        cvode_lib '/libsundials_cvode.a ' cvode_lib '/libsundials_fsunnonlinsolfixedpoint_mod.a"'];
else
    CVODE_include = '';
    CVODE = '';
    BUILD = [BUILD '_no_CVODE'];
    BUILD_MagTenseMicroMag = [BUILD_MagTenseMicroMag '_no_CVODE'];
end

if (VS_STUDIO)
    BUILD = [BUILD ' '];
    BUILD_MagTenseMicroMag = [BUILD_MagTenseMicroMag ' '];
else
    BUILD = ' ';
    BUILD_MagTenseMicroMag = ' ';
end


if (ispc)
    DEFINES = '-R2018a';
    FFLAGS = 'COMPFLAGS="$COMPFLAGS /free /O3 /fpp /real-size:64 /Qopenmp /assume:nocc_omp /fpe:0 /fp:source /nologo"';
    if (USE_CUDA)
        FFLAGS = [FFLAGS(1:(end-1)) ' /libs:static"'];
    end
    INCLUDE = ['-I' mkl_include ' -I' mkl_lp64 ' -I' NumericalIntegration_path BUILD '-I' DemagField_path BUILD '-I' TileDemagTensor_path ...
        BUILD '-I' MagTenseMicroMag_path BUILD_MagTenseMicroMag '-I' ForceIntegrator_path BUILD CVODE_include];
    LIBS = ['-L' MagTenseMicroMag_path BUILD_MagTenseMicroMag '-lMagTenseMicroMag -L' DemagField_path BUILD ...
        ' -lDemagField -L' TileDemagTensor_path BUILD ' -lTileDemagTensor -L' NumericalIntegration_path BUILD ...
        ' -lNumericalIntegration -L' ForceIntegrator_path BUILD ' -lMagneticForceIntegrator'];
    
    if (MKL_STATIC)
        MKL = ['-L' mkl_lib ' -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -lmkl_blas95_lp64 -llibiomp5md'];
    else
        MKL = ['-L' mkl_lib ' -lmkl_rt -lmkl_blas95_lp64'];
    end
else
    DEFINES = ['FC="' compiler_root '/bin/ifx" DEFINES="-DMATLAB_DEFAULT_RELEASE=R2018a"'];
    INCLUDE = ['INCLUDE="$INCLUDE -I' mkl_root '/include -I' NumericalIntegration_path ' -I' DemagField_path ...
        ' -I' TileDemagTensor_path ' -I' MagTenseMicroMag_path ' -I' ForceIntegrator_path ' -I' mkl_root];
    LIBS = ['-L' MagTenseMicroMag_path ' -lMagTenseMicroMag -L' DemagField_path ' -lDemagField -L' ...
        TileDemagTensor_path ' -lTileDemagTensor -L' NumericalIntegration_path ' -lNumericalIntegration -L' ...
        ForceIntegrator_path ' -lMagneticForceIntegrator'];
    FFLAGS = 'FFLAGS="';
    if (MKL_STATIC)
        MKL = ['-liomp5 -lpthread -lm -ldl -lifcoremt LINKLIBS="$LINKLIBS ' mkl_lib '/libmkl_intel_thread.a ' mkl_lib ...
            '/libmkl_core.a ' mkl_lib '/libmkl_blas95_lp64.a ' mkl_lib '/libmkl_intel_lp64.a"'];
        INCLUDE = [INCLUDE '/include/intel64/lp64 '];
    else
        MKL = ['-L' mkl_lib ' -lmkl_rt -lpthread -lm -ldl '];
        if (USE_CUDA)
            INCLUDE = [INCLUDE '/include/intel64/lp64 '];
            MKL = [MKL '-lmkl_blas95_lp64'];
        else
            FFLAGS = [FFLAGS '-i8 '];
            INCLUDE = [INCLUDE '/include/intel64/ilp64 '];
            MKL = [MKL '-lmkl_blas95_ilp64'];
        end
    end
    INCLUDE = [INCLUDE CVODE_include '"'];
    FFLAGS = [FFLAGS '-O3 -fpp -real-size 64 -qopenmp -assume nocc_omp -fpe0 -fp-model=source -fPIC -nologo -diag-disable 10006"'];
end

%%------------------------------------------------------------------
%%--------------- Build the MEX files ------------------------------
%%----------------------------------- ------------------------------
names = ["MagTenseLandauLifshitzSolver", "IterateMagnetization", "getHFromTiles", "getNFromTile", "getMagForce"]; 
if (~USE_CUDA)
    names(1) = "MagTenseLandauLifshitzSolverNoCUDA";
end
 
for i = 1:length(names)
    if names(i) == "MagTenseLandauLifshitzSolverNoCUDA"
        source = [mex_root "MagTenseLandauLifshitzSolver_mex.f90"];
        orig_name = "MagTenseLandauLifshitzSolver";
    else
        source = [mex_root names(i) '_mex.f90'];
        orig_name = names(i);
    end
    mex_str = ['mex' DEBUG DEFINES FFLAGS INCLUDE OBJS LIBS MKL CUDA CVODE join(source, '')];

    disp(join(mex_str, ' '))
    eval_MEX(join(mex_str, ' '))
    pause(pause_time)
    movefile(join([orig_name '_mex.mex' mex_suffix '64'], ''), join(['MEX_files/' names(i) '_mex.mex' mex_suffix '64'], ''));

    %--- Move the debug pdb files as well
    if (~USE_RELEASE)
        pause(pause_time)
        movefile(join([orig_name '_mex.mex' mex_suffix '64.pdb'], ''), join(['MEX_files/' names(i) '_mex.mex' mex_suffix '64.pdb'], ''));
    end
end

end

function eval_MEX(mex_str)
try
    eval(mex_str);
catch ME
    if (length(ME.message(:)) > 117)
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
    else
        disp(ME.message)
        rethrow(ME)
    end
end
end
