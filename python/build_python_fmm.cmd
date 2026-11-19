@echo off
rem Build the MagTense Python extension on Windows with the dip-fmm
rem demagnetization backend (USE_CDFMM=1).
rem
rem   python\build_python_fmm.cmd [clean]
rem
rem   clean  run `make clean` first; needed whenever the USE_* flags differ from
rem          those of the previous build.
rem
rem The root Makefile builds and installs dip-fmm itself (its cdfmm target). The
rem dip-fmm checkout is %CDFMM_DIR%, by default dip-fmm inside this repository
rem when there is one, otherwise dip-fmm next to it. Its tools\windows\env.cmd
rem sets up the build environment: the conda environment %CDFMM_CONDA_ENV%
rem (default: magtense-env, or the active one when it holds the Intel compilers)
rem and Visual Studio 2022.
setlocal
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

rem make is called by its absolute path: recipes run in conda's msys sh, where
rem Library\bin is shadowed and a recursive `make` is not found.
set "MAKE_EXE=%CDFMM_ENV%\Library\bin\make.exe"
rem The POSIX spelling of Library\bin survives msys' PATH translation, so ifx
rem and icx are also found from inside a recipe shell.
set "ENV_POSIX=%CDFMM_ENV:\=/%"
set "ENV_POSIX=/%ENV_POSIX::=%"
set "PATH=%ENV_POSIX%/Library/bin;%PATH%"
set "MAKE_ARGS=USE_CUDA=1 USE_CVODE=0 USE_MATLAB=0 USE_CDFMM=1 CDFMM_DIR=%CDFMM_DIR%"

cd /d "%MT_DIR%"
if /i "%~1"=="clean" (
  echo === make clean
  "%MAKE_EXE%" clean
)
echo === make python %MAKE_ARGS%
"%MAKE_EXE%" python %MAKE_ARGS%
if errorlevel 1 (echo MAGTENSE_BUILD_FAILED & exit /b 1)
echo MAGTENSE_BUILD_OK
