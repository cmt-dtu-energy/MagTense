@echo off

:: Preparation of terminal to find ifort and icl
:: Assuming there is a setvars.bat for Intel OneAPI on Windows
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

:: Navigate to wrapper file for python module
cd %RECIPE_DIR%\..\python\src\magtense\lib
del /Q /F %RECIPE_DIR%\..\python\src\magtense.egg-info
del /Q /F %RECIPE_DIR%\..\python\build
make clean
make

:: Navigate to distribution script
cd %RECIPE_DIR%\..\python
%PYTHON% -m pip install . -vvv --use-deprecated=legacy-resolver --no-deps --ignore-installed --no-cache-dir
