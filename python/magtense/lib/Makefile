#=======================================================================
#                   define the compiler names
#=======================================================================
# F90 = gfortran
F90 = ifort

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

ifeq ($(F90),ifort)
	F90FLAGS = /nologo /O3 /assume:nocc_omp /Qopenmp /real-size:64 /fp:precise /libs:static /threads
    FCOMP = intelvem # for f2py
	EXTRA_OPT = /assume:nounderscore /names:uppercase
endif

VPATH = ../../../source/NumericalIntegration/NumericalIntegration:\
../../../source/TileDemagTensor/TileDemagTensor:\
../../../source/DemagField/DemagField:\
../../../source/MagTenseMicroMag

CC = gcc
CFLAGS = -fPIC # ==> universal for ifort, gfortran, pgi
	
ifeq ($(OS),Windows_NT)
	MKLROOT="C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
	MKL= -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 \
	-lmkl_blas95_lp64 -lmkl_intel_thread -lmkl_core

	MKLROOT_COMP="C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\compiler"
	MKL_OMP=-L${MKLROOT_COMP}/lib/intel64_win -llibiomp5md

	LIBTOOL = LIB
	CLEAN_FILES = -del libsrc.lib *.o *.mod *.pyd
	CLEAN_FOLDER = -rmdir /s /q Release 2>NUL || VER>NUL
	CLEAN_FOLDER_WIN = -for /f %%i in ('dir /a:d /b src.win-*') do rd /s /q %%i
else
	MKLROOT=/opt/intel/oneapi/mkl/latest
	MKL= -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 \
	-lmkl_blas95_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl

	MKL_OMP=-L/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/ -liomp5

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o
	else
  		LIBTOOL = ar src
	endif
	CLEAN_FILES = -rm libsrc.a *.o *.mod *.so 
	CLEAN_FOLDER = -rm -r src.*
	CLEAN_FOLDER_WIN = -true
endif

OPT= ${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
OPT_PATH1 = ${MKLROOT}/include/intel64/lp64
OPT_PATH2 = ${MKLROOT}/include

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

all: _${PYTHON_MODN}.dll

clean:
	${CLEAN_FILES}
	${CLEAN_FOLDER}
	${CLEAN_FOLDER_WIN}

.f90.o:
	${F90} ${F90FLAGS} -I $(OPT) -c $< -o $@

.c.o:
	${CC} ${CFLAGS} -c $< -o $@

libsrc.lib: ${LIBSRC_OBJECTS}
	${LIBTOOL} $? /OUT:$@ 

_${PYTHON_MODN}.dll: libsrc.lib
	f2py -DNO_APPEND_FORTRAN -DUPPERCASE_FORTRAN --no-lower --build-dir . --fcompiler=${FCOMP} --opt='${F90FLAGS}' --opt='$(EXTRA_OPT)' -c -m ${PYTHON_MODN} -I:${OPT_PATH1} -I:${OPT_PATH2} -L. -llibsrc ${LIBSRC_WRAP_FILES} ${LIBSRC_OBJECTS} $(MKL) $(MKL_OMP)
