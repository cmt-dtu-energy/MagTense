@echo off

:: Navigate to wrapper file for python module
cd %RECIPE_DIR%\..\python\src\magtense\lib
del *.lib
rmdir /S /Q build
rmdir /S /Q %RECIPE_DIR%\..\python\src\magtense.egg-info
rmdir /S /Q %RECIPE_DIR%\..\python\build

:: Navigate to distribution script
cd %RECIPE_DIR%\..\python
%PYTHON% -m pip install . -vvv --use-deprecated=legacy-resolver --no-deps --ignore-installed --no-cache-dir
