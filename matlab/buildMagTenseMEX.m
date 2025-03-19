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
    mkl_lp64 = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\mkl\intel64\lp64"';
    mkl_lib = '"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib"';
    cuda_root = '"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\lib\x64"';
    cvode_include = '"C:\Program Files (x86)\sundials-4.1.0\instdir\fortran"';
    cvode_lib = '"C:\Program Files (x86)\sundials-4.1.0\instdir\lib"';
    mex_suffix = 'w';
else
    VS_STUDIO = false;
    MKL_STATIC = true;
    
    %--- The CONDA installation directory
    pre_str = getenv('CONDA_PREFIX');

    if (isempty(pre_str))
        %--- Check if miniconda is in /usr/share or in the users home
        share_dir = system('ls /usr/share/miniconda/envs/magtense-env/');
    
        if (share_dir == 0)
            pre_str = '/usr/share/miniconda';
        else
            %--- Get the username of the current user, which is where the miniconda is installed
            [~,username] = system('whoami');
            user = username(1:(end-1));
            pre_str = ['/home/' user '/miniconda3'];
        end
    end

    compiler_root = [pre_str '/envs/magtense-env'];
    mkl_root = [pre_str '/envs/magtense-env'];
    mkl_lib = [pre_str '/envs/magtense-env/lib'];
    cuda_root = [pre_str '/envs/magtense-env/lib'];
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
    if (ispc)
        CUDA = ['-L' cuda_root ' -lcublas -lcudart -lcuda -lcusparse'];
        OBJS = ['OBJS="$OBJS ' FortranCuda_path '/MagTenseCudaBlasICLWrapper.obj ' FortranCuda_path '/MagTenseCudaBlas.obj" '];
    else
        CUDA = ['-Wl,-Bdynamic ' '-L' cuda_root ' -lcublas -lcudart -lcusparse'];
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
    CVODE = ['-L' cvode_lib ' -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense' ...
    ' -lsundials_fnvecserial_mod -lsundials_cvode -lsundials_fsunnonlinsolfixedpoint_mod'];
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
    FFLAGS = 'COMPFLAGS="$COMPFLAGS /free /nologo /real-size:64 /O2 /assume:nocc_omp /Qopenmp /fpp /fpe:0 /fp:source /fp:precise"';
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
    INCLUDE = ['INCLUDE="$INCLUDE -I' NumericalIntegration_path ' -I' DemagField_path ...
        ' -I' TileDemagTensor_path ' -I' MagTenseMicroMag_path ' -I' ForceIntegrator_path ... 
        ' -I' mkl_root '/include' ' -I' mkl_root '/include/intel64/lp64'];
    LIBS = ['LINKLIBS=''$LINKLIBS ' '-L' MagTenseMicroMag_path ' -lMagTenseMicroMag -L' ...
        ForceIntegrator_path ' -lMagneticForceIntegrator -L' DemagField_path ' -lDemagField -L' ...
        TileDemagTensor_path ' -lTileDemagTensor -L' NumericalIntegration_path ' -lNumericalIntegration'''];
    if (MKL_STATIC)
        LIBS = [LIBS(1:(end-1)) ' ' mkl_lib '/libmkl_blas95_lp64.a -Wl,--start-group ' ...
               mkl_lib '/libmkl_intel_lp64.a ' mkl_lib '/libmkl_intel_thread.a ' ...
               mkl_lib '/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -static-intel'''];
        MKL = [];
        %INCLUDE = [INCLUDE '/include/intel64/lp64 '];
        if (USE_CUDA)
            LIBS = [LIBS(1:(end-1)) ' ' CUDA ''''];
            CUDA = '';
        end
    else
        MKL = ['-L' mkl_root '/lib/intel64 -lmkl_rt -lpthread -lm -ldl '];
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
    FFLAGS = ['FFLAGS="-O3 -fpp -real-size 64 -fpe0 -fp-model=source -fPIC -nologo -diag-disable 10006"'];
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