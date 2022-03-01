#=======================================================================
#                   define the compiler names
#=======================================================================
F90 = gfortran

#=======================================================================
#                     additional flags
#=======================================================================
ifeq ($(F90),gfortran)
	F90FLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512
	ifeq ($(OS),Windows_NT)
    	FCOMP = gnu95 --compiler=mingw32
	else
		FCOMP = gnu95
	endif
endif

MKLROOT=/opt/intel/oneapi/mkl/latest
OPT= ${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
MKL= -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 \
-lmkl_blas95_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl
MKL_OMP=-L/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/ -liomp5

VPATH = ../../../source/NumericalIntegration/NumericalIntegration:\
../../../source/TileDemagTensor/TileDemagTensor:\
../../../source/DemagField/DemagField:\
../../../source/MagTenseMicroMag

FPP = gfortran -E # gfortran f90wrap temp files only, not compilation
FPP_F90FLAGS = -x f95-cpp-input -fPIC
CC = gcc
CFLAGS = -fPIC # ==> universal for ifort, gfortran, pgi

ifeq ($(OS),Windows_NT)
	LIBTOOL = ar src
	CLEAN_FILES = -del libsrc.a *.o *.mod *.pyd
	CLEAN_FOLDER = -rmdir /s /q Release 2>NUL || VER>NUL
	CLEAN_FOLDER_WIN = -for /f %%i in ('dir /a:d /b src.win-*') do rd /s /q %%i
else
	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o
	else
  		LIBTOOL = ar src
	endif
	CLEAN_FILES = -rm libsrc.a *.o *.mod *.so 
	CLEAN_FOLDER = -rm -r src.*
	CLEAN_FOLDER_WIN = -true
endif

# ======================================================================
# PROJECT CONFIG, do not put spaced behind the variables
# ======================================================================
# Python module name
PYTHON_MODN = magtensesource

#=======================================================================
#       List all source files required for the project
#=======================================================================
# names (without suffix), f90 sources
LIBSRC_SOURCES = IntegrationDataTypes quadpack SpecialFunctions TileTensorHelperFunctions \
TileRectangularPrismTensor TileCircPieceTensor TileCylPieceTensor TilePlanarCoilTensor TileTriangle \
TileNComponents DemagFieldGetSolution spline MagParameters IterateMagnetSolution \
MicroMagParameters RKSuitef90 ODE_Solvers LLODE_Debug util MagTenseMicroMagPyIO LandauLifshitzEquationSolver

# file names
LIBSRC_FILES = $(addsuffix .f90,${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(addsuffix .o,${LIBSRC_SOURCES})

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================
# names (without suffix), f90 sources
LIBSRC_WRAP_SOURCES = FortranToPythonIO

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

#=======================================================================
#                 Relevant suffixes
#=======================================================================
.SUFFIXES: .f90

#=======================================================================
#
#=======================================================================
.PHONY: all clean

all: _${PYTHON_MODN}.so

clean:
	${CLEAN_FILES}
	${CLEAN_FOLDER}
	${CLEAN_FOLDER_WIN}

.f90.o:
	${F90} ${F90FLAGS} -I $(OPT) -c $< -o $@

.c.o:
	${CC} ${CFLAGS} -c $< -o $@

libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?

_${PYTHON_MODN}.so: libsrc.a
	f2py --build-dir . --fcompiler=${FCOMP} --opt='${F90FLAGS}' -c -m ${PYTHON_MODN} -L. -lsrc ${LIBSRC_WRAP_FILES} ${LIBSRC_OBJECTS} $(MKL) $(MKL_OMP)
