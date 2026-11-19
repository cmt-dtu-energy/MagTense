@echo off
rem Build the MagTense MEX files on Windows with the dip-fmm demagnetization
rem backend, with the root Makefile (not the Visual Studio solution).
rem
rem   matlab\build_mex_fmm.cmd [cuda|nocuda] [mexonly]      (default: cuda)
rem
rem Step 1 (conda environment): make clean, then the libraries with USE_MATLAB=1
rem   and USE_CDFMM=1. dip-fmm is rebuilt with static oneMKL, so that
rem   cdfmm_c.dll does not bind to the MKL inside MATLAB, and always with its
rem   CUDA backend (CDFMM_ENABLE_CUDA=ON), so that the one cdfmm_c.dll in
rem   MEX_files serves both the CUDA and the NoCUDA solver MEX.
rem Step 2 (plain PATH + standalone oneAPI ifx): buildMagTenseMEX in MATLAB.
rem   `mexonly` runs only this step, on libraries that are already built.
rem
rem CDFMM_DIR        the dip-fmm checkout; default: dip-fmm inside this repository
rem                  when there is one, otherwise dip-fmm next to it
rem CDFMM_CONDA_ENV  the conda environment; default: magtense-env
rem MATLAB_ROOT      default: C:\Program Files\MATLAB\R2024b
rem IFX_BIN          folder of the ifx that MATLAB's mex is set up with; default:
rem                  C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin
rem
rem NOTE: make clean removes the Python extension; run python\build_python_fmm.cmd
rem clean afterwards to get it back.
setlocal
set "VARIANT=%~1"
if "%VARIANT%"=="" set "VARIANT=cuda"
set "USE_CUDA=1"
set "ML_CUDA=true"
if /i "%VARIANT%"=="nocuda" set "USE_CUDA=0"
if /i "%VARIANT%"=="nocuda" set "ML_CUDA=false"

set "ORIG_PATH=%PATH%"
for %%i in ("%~dp0..") do set "MT_DIR=%%~fi"
if not defined CDFMM_DIR if exist "%MT_DIR%\dip-fmm\CMakeLists.txt" set "CDFMM_DIR=%MT_DIR%\dip-fmm"
if not defined CDFMM_DIR for %%i in ("%MT_DIR%\..\dip-fmm") do set "CDFMM_DIR=%%~fi"
if not exist "%CDFMM_DIR%\tools\windows\env.cmd" (
  echo No dip-fmm checkout with tools\windows\env.cmd at "%CDFMM_DIR%"
  exit /b 1
)
if not defined CDFMM_CONDA_ENV set "CDFMM_CONDA_ENV=magtense-env"
call "%CDFMM_DIR%\tools\windows\env.cmd"
if errorlevel 1 exit /b 1
@echo off
set "CDFMM_DIR=%CDFMM_DIR:\=/%"
if not defined MATLAB_ROOT set "MATLAB_ROOT=C:\Program Files\MATLAB\R2024b"
if not defined IFX_BIN set "IFX_BIN=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"

set "MAKE_EXE=%CDFMM_ENV%\Library\bin\make.exe"
set "ENV_POSIX=%CDFMM_ENV:\=/%"
set "ENV_POSIX=/%ENV_POSIX::=%"
set "PATH=%ENV_POSIX%/Library/bin;%PATH%"

if /i "%~2"=="mexonly" goto mexstep

cd /d "%MT_DIR%"
echo === make clean
"%MAKE_EXE%" clean
echo === make all (USE_MATLAB=1 USE_CDFMM=1 USE_CUDA=%USE_CUDA%)
"%MAKE_EXE%" all USE_CUDA=%USE_CUDA% USE_CVODE=0 USE_MATLAB=1 "MATLAB_INCLUDE=%MATLAB_ROOT%\extern\include" USE_CDFMM=1 CDFMM_DIR=%CDFMM_DIR% CDFMM_ENABLE_CUDA=ON
if errorlevel 1 (echo MEX_LIBS_FAILED & exit /b 1)

:mexstep
echo === buildMagTenseMEX
rem MATLAB's mex needs the standalone oneAPI ifx on a Windows-form PATH.
set "PATH=%IFX_BIN%;%ORIG_PATH%"
cd /d "%MT_DIR%\matlab"
"%MATLAB_ROOT%\bin\matlab.exe" -batch "buildMagTenseMEX('USE_CUDA', %ML_CUDA%, 'USE_CVODE', false, 'USE_FMM', true, 'fmm_root', '%CDFMM_DIR%/local')"
if errorlevel 1 (echo MEX_BUILD_FAILED & exit /b 1)
echo MEX_BUILD_OK
