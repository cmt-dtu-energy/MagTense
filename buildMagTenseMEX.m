function buildMagTenseMEX(USE_RELEASE,USE_CUDA,USE_CVODE)
%use clear all as this also clears dependencies to the .mex files and thus they can be overwritten
      
clearvars -except USE_RELEASE USE_CUDA USE_CVODE
  
if ~exist('USE_RELEASE','var')
    USE_RELEASE = true;
end
if ~exist('USE_CUDA','var')
    USE_CUDA = true;
end
if ~exist('USE_CVODE','var')
    USE_CVODE = true;
end

%% The different possible cases
if (USE_RELEASE)
    Debug_flag = '' ;
    build_str = 'release';
else
    Debug_flag = '-g' ;
    build_str = 'debug';
end

if (USE_CUDA)
    CUDA_str  = '''-Lc:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.2\lib\x64\'' -lcublas -lcudart -lcuda -lcusparse';
    if (USE_RELEASE)
        build_str_MagTenseMicroMag = 'release';
    else
        build_str_MagTenseMicroMag = 'debug';
    end
else
    CUDA_str  = '';
    if (USE_RELEASE)
        build_str_MagTenseMicroMag = 'Release_no_CUDA';
    else
        build_str_MagTenseMicroMag = 'Debug_no_CUDA';
    end
end

if (USE_CVODE)
    CVODE_str = ['''-Lc:\Program Files (x86)\sundials-4.1.0\instdir\lib'' -lsundials_nvecserial ' ...
             '''-Lc:\Program Files (x86)\sundials-4.1.0\instdir\lib'' -lsundials_sunmatrixdense ' ...
             '''-Lc:\Program Files (x86)\sundials-4.1.0\instdir\lib'' -lsundials_sunlinsoldense ' ...
             '''-Lc:\Program Files (x86)\sundials-4.1.0\instdir\lib'' -lsundials_fnvecserial_mod ' ...
             '''-Lc:\Program Files (x86)\sundials-4.1.0\instdir\lib'' -lsundials_cvode ' ...
             '''-Lc:\Program Files (x86)\sundials-4.1.0\instdir\lib'' -lsundials_fsunnonlinsolfixedpoint_mod ' ...
             '''-Ic:\Program Files (x86)\sundials-4.1.0\instdir\include'''];
    if (USE_RELEASE)
        build_str_NumericalIntegration = 'release';
    else
        build_str_NumericalIntegration = 'debug';
    end
else
    CVODE_str  = '';
    if (USE_RELEASE)
        build_str_NumericalIntegration = 'Release_no_CVODE';
    else
        build_str_NumericalIntegration = 'Debug_no_CVODE';
    end
end

%--- The individual Fortran parts
MagneticForceIntegrator_str = ['-Lsource\MagneticForceIntegrator\MagneticForceIntegrator\x64\' build_str '\ -lMagneticForceIntegrator -Isource\MagneticForceIntegrator\MagneticForceIntegrator\x64\' build_str '\'];
DemagField_str              = ['-Lsource\DemagField\DemagField\x64\' build_str '\ -lDemagField -Isource\DemagField\DemagField\x64\' build_str '\'];
NumericalIntegration_str    = ['-Lsource\NumericalIntegration\NumericalIntegration\x64\' build_str_NumericalIntegration '\ -lNumericalIntegration -Isource\NumericalIntegration\NumericalIntegration\x64\' build_str_NumericalIntegration '\'];
TileDemagTensor_str         = ['-Lsource\TileDemagTensor\TileDemagTensor\x64\' build_str '\ -lTileDemagTensor -Isource\TileDemagTensor\TileDemagTensor\x64\' build_str '\'];
MagTenseMicroMag_str        = ['-Lsource\MagTenseMicroMag\x64\' build_str_MagTenseMicroMag '\ -lMagTenseMicroMag -Isource\MagTenseMicroMag\x64\' build_str_MagTenseMicroMag '\'];

Options_str                 = 'COMPFLAGS="$COMPFLAGS /O3 /extend_source:132 /real_size:64 /fpe:0" -R2018a';
MKL_str                     = '''-Lc:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019\windows\mkl\lib\intel64_win\'' -lmkl_intel_lp64 -lmkl_blas95_lp64 ''-Ic:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019\windows\mkl\include\''';

%%--------------------------------------------------------------------------------------------------------
%%----------------------------------------------------- Build the MEX files ------------------------------
%%------------------------------------------------------------------------- ------------------------------
%% MagTenseLandauLifshitzSolver_mex
Source_str = 'source\MagTenseMEX\MagTenseMEX\MagTenseLandauLifshitzSolver_mex.f90';
mex_str = ['mex ' Debug_flag ' ' NumericalIntegration_str ' ' DemagField_str ' ' TileDemagTensor_str ' ' CUDA_str ' ' CVODE_str ' ' MKL_str ' ' MagTenseMicroMag_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    movefile MagTenseLandauLifshitzSolver_mex.mexw64.pdb matlab\MEX_files\MagTenseLandauLifshitzSolver_mex.mexw64.pdb
end
movefile MagTenseLandauLifshitzSolver_mex.mexw64 matlab\MEX_files\MagTenseLandauLifshitzSolver_mex.mexw64

%% IterateMagnetization_mex
Source_str = 'source\MagTenseMEX\MagTenseMEX\IterateMagnetization_mex.f90';
mex_str = ['mex ' Debug_flag ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    movefile IterateMagnetization_mex.mexw64.pdb matlab\MEX_files\IterateMagnetization_mex.mexw64.pdb
end
movefile IterateMagnetization_mex.mexw64 matlab\MEX_files\IterateMagnetization_mex.mexw64
        
%% getHFromTiles_mex
Source_str = 'source\MagTenseMEX\MagTenseMEX\getHFromTiles_mex.f90';
mex_str = ['mex ' Debug_flag ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    movefile getHFromTiles_mex.mexw64.pdb matlab\MEX_files\getHFromTiles_mex.mexw64.pdb
end
movefile getHFromTiles_mex.mexw64 matlab\MEX_files\getHFromTiles_mex.mexw64
    
%% getNFromTile_mex
Source_str = 'source\MagTenseMEX\MagTenseMEX\getNFromTile_mex.f90';
mex_str = ['mex ' Debug_flag ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    movefile getNFromTile_mex.mexw64.pdb matlab\MEX_files\getNFromTile_mex.mexw64.pdb
end
movefile getNFromTile_mex.mexw64 matlab\MEX_files\getNFromTile_mex.mexw64
    
%% getMagForce_mex
Source_str = 'source\MagTenseMEX\MagTenseMEX\getMagForce_mex.f90';
mex_str = ['mex ' Debug_flag ' ' MagneticForceIntegrator_str ' ' DemagField_str ' ' NumericalIntegration_str ' ' TileDemagTensor_str ' ' Source_str ' ' Options_str];
eval(mex_str) 
if ~USE_RELEASE
    movefile getMagForce_mex.mexw64.pdb matlab\MEX_files\getMagForce_mex.mexw64.pdb
end
movefile getMagForce_mex.mexw64 matlab\MEX_files\getMagForce_mex.mexw64
    
end
