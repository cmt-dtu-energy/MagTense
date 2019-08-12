# Creating an importable Python module from Fortran source code (**.pyd** on Windows / **.so** on Linux)

Requirements:

- Python installation
- Numpy package
- Fortran compiler: gfortran or ifort

## Fortran compiler: **gfortran** from MinGW64

Open terminal, navigate to this folder and run Makefile in activated Python virtual environment.

- Windows:

    ```Powershell
    cd PATH\TO\THIS\FOLDER
    ..\venv\Scripts\activate.bat
    mingw32-make.exe
    ```

- Linux

    ```bash
    cd PATH/TO/THIS/FOLDER
    source ../venv/bin/activate
    make
    ```

## Fortran file with functions to wrap for Python

- FortranToPythonIO.f90

## Fortran compiler: **ifort** compiler from Intel Parallel Studio XE - CURRENTLY NOT WORKING

Change compiler in Makefile to ifort:

```Makefile
#F90 = gfortran
F90 = ifort
```

Open terminal, activate Intel Software Environment

```Powershell
cd C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019\windows\bin\
CALL ifortvars.bat intel64
```

Navigate to this folder, activate Python virtual environment and run Makefile.

-----> LINKING ERRORS
