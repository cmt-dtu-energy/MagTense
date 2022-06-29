function buildMagTenseMEX(USE_RELEASE,USE_CUDA,USE_CVODE)
%use clear all as this also clears dependencies to the .mex files and thus they can be overwritten
      
arguments
    USE_RELEASE {mustBeNumericOrLogical} = true;
    USE_CUDA {mustBeNumericOrLogical} = true;
    USE_CVODE {mustBeNumericOrLogical} = true;
end

pause_time = 1; %Time to wait between making and moving the generated files

%% The different possible cases
if (ispc)
    compiler_str = '';
    MagneticForceIntegrator_lib_str = '-lMagneticForceIntegrator';
    DemagField_lib_str              = '-lDemagField';
    NumericalIntegration_lib_str    = '-lNumericalIntegration';
    TileDemagTensor_lib_str         = '-lTileDemagTensor';
    MagTenseMicroMag_lib_str        = '-lMagTenseMicroMag';
    MEX_str = 'w';
    % Necessary for Matlab to recognize oneAPI
    setenv('IFORT_COMPILER19','C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows');
else
    %compiler_str = 'FC=''/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort -fpp -fpic -free -Wl,--no-as-needed -mkl -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -static-intel ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_blas95_lp64.a'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_lapack95_lp64.a'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_lp64.a'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_thread.a'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_core.a'' -liomp5 -lpthread -lm -ldl ''-I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include/intel64/lp64/'' ''-I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include''''';
    

    compiler_str = 'FC=''/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort -fopenmp -fpp -fpic -free -mkl -static-intel -I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include/intel64/lp64 -I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include'''; 
    %compiler_str = 'FC=''/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort -fpp -fpic -free -mkl /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_blas95_lp64.a  -liomp5 -lpthread -lm -ldl''';    
    %compiler_str = 'FC=''/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort -fpp -fpic -free -mkl -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -static-intel -libs:static -I../NumericalIntegration/NumericalIntegration/ -I../DemagField/DemagField/ -I/apps/external/matlab/2019/b/extern/include/ -I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/ -I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include/ -I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include/intel64/lp64 -I../TileDemagTensor/TileDemagTensor/ -Wl, --start-group /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_lp64.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_blas95_lp64.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_lapack95_lp64.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_thread.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl''';
    %compiler_str = 'FC=''/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort -fpp -fpic -free -mkl -static-intel Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64  -Wl, --start-group /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_lp64.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_blas95_lp64.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_lapack95_lp64.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_thread.a /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl''';
    MagneticForceIntegrator_lib_str = '-lMagneticForceIntegrator ''-Lsource/MagneticForceIntegrator/MagneticForceIntegrator/MagneticForceIntegratorLib.a''';
    DemagField_lib_str              = '-lDemagField ''-Lsource/DemagField/DemagField/DemagFieldLib.a''';
    NumericalIntegration_lib_str    = '-lNumericalIntegration ''-Lsource/NumericalIntegration/NumericalIntegration/NumericalIntegrationLib.a''';
    TileDemagTensor_lib_str         = '-lTileDemagTensor ''-Lsource/TileDemagTensor/TileDemagTensor/TileDemagTensorLib.a''';
    MagTenseMicroMag_lib_str        = '-lMagTenseMicroMag ''-Lsource/MagTenseMicroMag/MagTenseMicroMagLib.a''';
    MEX_str = 'a';
end

if (ispc)
    if (USE_RELEASE)
        Debug_flag = '' ;
        build_str = 'x64/Release';
    else
        Debug_flag = '-g' ;
        build_str = 'x64/Debug';
    end
else
    if (USE_RELEASE)
        Debug_flag = '-g' ;
    else
        Debug_flag = '-g' ;
    end 
    build_str = '';
end

if (ispc)
    if (USE_CUDA)
        CUDA_str  = '''-Lc:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.6/lib/x64/'' -lcublas -lcudart -lcuda -lcusparse';
        build_str_MagTenseMicroMag = build_str;
        build_str_NO_CUDA_MagTenseMicroMag = [build_str '_no_CUDA'];
    else
        CUDA_str  = '';        
        build_str_MagTenseMicroMag = [build_str '_no_CUDA'];
    end
else
    if (USE_CUDA)
        CUDA_str  = '';
    else
        CUDA_str  = '';
    end
    build_str_MagTenseMicroMag = '';
end

if (ispc)
    if (USE_CVODE)
        CVODE_str = ['''-Lc:/Program Files (x86)/sundials-4.1.0/instdir/lib'' -lsundials_nvecserial ' ...
                 '''-Lc:/Program Files (x86)/sundials-4.1.0/instdir/lib'' -lsundials_sunmatrixdense ' ...
                 '''-Lc:/Program Files (x86)/sundials-4.1.0/instdir/lib'' -lsundials_sunlinsoldense ' ...
                 '''-Lc:/Program Files (x86)/sundials-4.1.0/instdir/lib'' -lsundials_fnvecserial_mod ' ...
                 '''-Lc:/Program Files (x86)/sundials-4.1.0/instdir/lib'' -lsundials_cvode ' ...
                 '''-Lc:/Program Files (x86)/sundials-4.1.0/instdir/lib'' -lsundials_fsunnonlinsolfixedpoint_mod ' ...
                 '''-Ic:/Program Files (x86)/sundials-4.1.0/instdir/include'''];
    else
        CVODE_str  = '';
        build_str = [build_str '_no_CVODE'];
        build_str_MagTenseMicroMag = [build_str_MagTenseMicroMag '_no_CVODE'];
    end
else
    if (USE_CVODE)
        CVODE_str = [''];
    else
        CVODE_str  = '';
    end
end

%--- The individual Fortran parts
MagneticForceIntegrator_str = ['-Lsource/MagneticForceIntegrator/MagneticForceIntegrator/' build_str '/ ' MagneticForceIntegrator_lib_str ' -Isource/MagneticForceIntegrator/MagneticForceIntegrator/' build_str '/'];
DemagField_str              = ['-Lsource/DemagField/DemagField/' build_str '/ ' DemagField_lib_str ' -Isource/DemagField/DemagField/' build_str '/'];
NumericalIntegration_str    = ['-Lsource/NumericalIntegration/NumericalIntegration/' build_str '/ ' NumericalIntegration_lib_str ' -Isource/NumericalIntegration/NumericalIntegration/' build_str '/'];
TileDemagTensor_str         = ['-Lsource/TileDemagTensor/TileDemagTensor/' build_str '/ ' TileDemagTensor_lib_str ' -Isource/TileDemagTensor/TileDemagTensor/' build_str '/'];
MagTenseMicroMag_str        = ['-Lsource/MagTenseMicroMag/' build_str_MagTenseMicroMag '/ ' MagTenseMicroMag_lib_str ' -Isource/MagTenseMicroMag/' build_str_MagTenseMicroMag '/'];

Options_str                 = 'COMPFLAGS="$COMPFLAGS /Qopenmp /O2 /free /real-size:64 /fpe:0" -R2018a'; %/real_size:64 %/O3
if (ispc)
    MKL_str                 = '''-LC:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/intel64'' -lmkl_intel_lp64 -lmkl_blas95_lp64 ''-IC:\Program Files (x86)\Intel\oneAPI\mkl\latest\include''';%''-Ic:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019/windows/mkl/include/''';
else
    %MKL_str                 = '''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/'' -lmkl_intel_lp64 -lmkl_blas95_lp64 ''-I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include/''';
    MKL_str                 = '/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_blas95_lp64.a  -liomp5 -lpthread -lm -ldl';%'''-I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/'' ''-I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include/'' ''-I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/include/intel64/lp64''  ''-I/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_blas95_lp64.a'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_lp64.a'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_intel_thread.a'' ''-L/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64/libmkl_core.a''';
end

%%--------------------------------------------------------------------------------------------------------
%%----------------------------------------------------- Build the MEX files ------------------------------
%%------------------------------------------------------------------------- ------------------------------
%% MagTenseLandauLifshitzSolver_mex
Source_str = 'source/MagTenseMEX/MagTenseMEX/MagTenseLandauLifshitzSolver_mex.f90';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' MagTenseMicroMag_str ' ' DemagField_str ' ' TileDemagTensor_str ' ' NumericalIntegration_str ' ' CUDA_str ' ' CVODE_str ' ' MKL_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    pause(pause_time)
    movefile(['MagTenseLandauLifshitzSolver_mex.mex' MEX_str '64.pdb'],['matlab/MEX_files/MagTenseLandauLifshitzSolver_mex.mex' MEX_str '64.pdb']);
end
pause(pause_time)
movefile(['MagTenseLandauLifshitzSolver_mex.mex' MEX_str '64'],['matlab/MEX_files/MagTenseLandauLifshitzSolver_mex.mex' MEX_str '64']);

%% No CUDA version of MagTenseLandauLifshitzSolver_mex
MagTenseMicroMag_str        = ['-Lsource/MagTenseMicroMag/' build_str_NO_CUDA_MagTenseMicroMag '/ ' MagTenseMicroMag_lib_str ' -Isource/MagTenseMicroMag/' build_str_MagTenseMicroMag '/'];
Source_str = 'source/MagTenseMEX/MagTenseMEX/MagTenseLandauLifshitzSolver_mex.f90';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' MagTenseMicroMag_str ' ' DemagField_str ' ' TileDemagTensor_str ' ' NumericalIntegration_str ' ' CVODE_str ' ' MKL_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    pause(pause_time)
    movefile(['MagTenseLandauLifshitzSolver_mex.mex' MEX_str '64.pdb'],['matlab/MEX_files/MagTenseLandauLifshitzSolverNoCUDA_mex.mex' MEX_str '64.pdb']);
end
pause(pause_time)
movefile(['MagTenseLandauLifshitzSolver_mex.mex' MEX_str '64'],['matlab/MEX_files/MagTenseLandauLifshitzSolverNoCUDA_mex.mex' MEX_str '64']);


%% IterateMagnetization_mex
Source_str = 'source/MagTenseMEX/MagTenseMEX/IterateMagnetization_mex.f90';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' DemagField_str ' ' TileDemagTensor_str ' ' NumericalIntegration_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    pause(pause_time)
    movefile(['IterateMagnetization_mex.mex' MEX_str '64.pdb'],['matlab/MEX_files/IterateMagnetization_mex.mex' MEX_str '64.pdb']);
end
pause(pause_time)
movefile(['IterateMagnetization_mex.mex' MEX_str '64'],['matlab/MEX_files/IterateMagnetization_mex.mex' MEX_str '64']);
        
%% getHFromTiles_mex
Source_str = 'source/MagTenseMEX/MagTenseMEX/getHFromTiles_mex.f90';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' DemagField_str ' ' TileDemagTensor_str ' ' NumericalIntegration_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    pause(pause_time)
    movefile(['getHFromTiles_mex.mex' MEX_str '64.pdb'],['matlab/MEX_files/getHFromTiles_mex.mex' MEX_str '64.pdb']);
end
pause(pause_time)
movefile(['getHFromTiles_mex.mex' MEX_str '64'],['matlab/MEX_files/getHFromTiles_mex.mex' MEX_str '64']);
    
%% getNFromTile_mex
Source_str = 'source/MagTenseMEX/MagTenseMEX/getNFromTile_mex.f90';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    pause(pause_time)
    movefile(['getNFromTile_mex.mex' MEX_str '64.pdb'],['matlab/MEX_files/getNFromTile_mex.mex' MEX_str '64.pdb']);
end
pause(pause_time)
movefile(['getNFromTile_mex.mex' MEX_str '64'],['matlab/MEX_files/getNFromTile_mex.mex' MEX_str '64']);
    
%% getMagForce_mex
Source_str = 'source/MagTenseMEX/MagTenseMEX/getMagForce_mex.f90';
mex_str = ['mex ' compiler_str ' ' Debug_flag ' ' MagneticForceIntegrator_str ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    pause(pause_time)
    movefile(['getMagForce_mex.mex' MEX_str '64.pdb'],['matlab/MEX_files/getMagForce_mex.mex' MEX_str '64.pdb']);
end
pause(pause_time)
movefile(['getMagForce_mex.mex' MEX_str '64'],['matlab/MEX_files/getMagForce_mex.mex' MEX_str '64']);
    
end
