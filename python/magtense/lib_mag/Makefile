#=======================================================================
#                   define the compiler names
#=======================================================================
F90 = gfortran
#F90 = ifort

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
endif

VPATH = ../../source/NumericalIntegration/NumericalIntegration:\
../../source/TileDemagTensor/TileDemagTensor:\
../../source/DemagField/DemagField
FPP = gfortran -E # gfortran f90wrap temp files only, not compilation
FPP_F90FLAGS = -x f95-cpp-input -fPIC
CC = gcc
CFLAGS = -fPIC # ==> universal for ifort, gfortran, pgi

ifeq ($(OS),Windows_NT)
	LIBTOOL = ar src
	CLEAN_FILES = -del libsrc.a *.o *.mod *.pyd
	CLEAN_FOLDER = -rmdir Release src.win-amd64-3.7 /s /Q
else
	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o
	else
  		LIBTOOL = ar src
	endif
	CLEAN_FILES = -rm libsrc.a *.o *.mod *.so 
	CLEAN_FOLDER = -rm -r src.*
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
TileRectangularPrismTensor TileCircPieceTensor TileCylPieceTensor TilePlanarCoilTensor TileTriangle\
TileNComponents DemagFieldGetSolution spline MagParameters IterateMagnetSolution

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

.f90.o:
	${F90} ${F90FLAGS} -c $< -o $@

.c.o:
	${CC} ${CFLAGS} -c $< -o $@

libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?

_${PYTHON_MODN}.so: libsrc.a
	f2py --build-dir . --fcompiler=${FCOMP} --opt='${F90FLAGS}' -c -m ${PYTHON_MODN} -L. -lsrc ${LIBSRC_WRAP_FILES} ${LIBSRC_OBJECTS}
