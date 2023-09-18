@echo off

:: Preparation of terminal to find ifort
:: Assuming there is a setvars.bat for Intel OneAPI on Windows
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

:: Navigate to wrapper file for python module
cd %RECIPE_DIR%\..\python\src\magtense\lib
make

:: Navigate to distribution script
cd %RECIPE_DIR%\..\python
%PYTHON% -m pip install . -vvv --use-deprecated=legacy-resolver --no-deps --ignore-installed --no-cache-dir

:: Clean up
cd %RECIPE_DIR%\..\python\src\magtense\lib
make clean
del /Q /F %RECIPE_DIR%\..\python\src\magtense.egg-info
del /Q /F %RECIPE_DIR%\..\python\build
