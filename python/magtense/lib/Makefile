#=======================================================================
#                       compiler names and flags
#=======================================================================
F90 = ifort
CC = gcc
CFLAGS = -fPIC
PYTHON_MODN = magtensesource
USE_CUDA = 1
USE_CVODE = 0

VPATH = ../../../source/NumericalIntegration/NumericalIntegration:\
../../../source/TileDemagTensor/TileDemagTensor:\
../../../source/DemagField/DemagField:\
../../../source/MagTenseMicroMag:\
../../../source/MagTenseFortranCuda/cuda


ifeq ($(OS),Windows_NT)
	F90FLAGS = /nologo /O3 /assume:nocc_omp /Qopenmp /real-size:64 /fp:precise \
	/threads /fpp /DUSE_CUDA=${USE_CUDA} /DUSE_CVODE=${USE_CVODE} /DUSE_MATLAB=0

	EXTRA_OPT = -DNO_APPEND_FORTRAN -DUPPERCASE_FORTRAN --no-lower --opt='${F90FLAGS} /assume:nounderscore /names:uppercase'
	FCOMP = intelvem

	ifeq ($(USE_CUDA),1)
		CUDA = -L${CONDA_PREFIX}/Lib/x64 -lcublas -lcudart -lcuda -lcusparse
		CUDA_OBJECTS = MagTenseCudaBlas.obj MagTenseCudaBlasICLWrapper.obj
	endif

	MKL = -L${CONDA_PREFIX}/Library/lib -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
	MKL_OMP = -L"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\compiler\lib\intel64_win" -llibiomp5md
	OPT = ${CONDA_PREFIX}/Library/include -I${CONDA_PREFIX}/Library/include/intel64/ilp64
	LIB_OPT = -I:${OPT} -L. -llibsrc
	
	LIBTOOL = LIB $? /OUT:$@
	STATIC_LIB = libsrc.lib
	CLEAN_FILES = -del libsrc.lib *.o *.mod *.pyd
	CLEAN_FOLDER = -rmdir /s /q Release 2>NUL || VER>NUL
	CLEAN_FOLDER_WIN = -for /f %%i in ('dir /a:d /b src.win-*') do rd /s /q %%i
	PYTHON_MODN_ALL = _${PYTHON_MODN}.dll

else
	export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib
	
	F90FLAGS = -fPIC -nologo -O3 -assume nocc_omp -qopenmp -real-size 64 \
	-fp-model precise -threads -fpp -DUSE_CUDA=${USE_CUDA} -DUSE_CVODE=${USE_CVODE} -DUSE_MATLAB=0

	EXTRA_OPT = --opt='${F90FLAGS}'
	FCOMP = intelem

	ifeq ($(USE_CUDA),1)
		CUDA = -L${CONDA_PREFIX}/lib -lcublas -lcudart -lcuda -lcusparse
		CUDA_OBJECTS = MagTenseCudaBlas.o MagTenseCudaBlasICLWrapper.o
	endif

	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -lpthread -lm -ldl -liomp5 -lmkl_blas95_lp64
	OPT = ${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/intel64/ilp64
	LIB_OPT = -I${OPT} -L. -lsrc

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	
	STATIC_LIB = libsrc.a
	CLEAN_FILES = -rm libsrc.a *.o *.mod *.so 
	CLEAN_FOLDER = -rm -r src.*
	CLEAN_FOLDER_WIN = -true
	PYTHON_MODN_ALL = _${PYTHON_MODN}.so
endif

ifeq ($(USE_CUDA),0)
	CUDA =
	CUDA_OBJECTS =
endif

#=======================================================================
#       List all source files required for the project
#=======================================================================
# names (without suffix), f90 sources
LIBSRC_SOURCES = IntegrationDataTypes quadpack SpecialFunctions TileTensorHelperFunctions \
TileRectangularPrismTensor TileCircPieceTensor TileCylPieceTensor TilePlanarCoilTensor \
TileTriangle TileNComponents DemagFieldGetSolution spline MagParameters IterateMagnetSolution \
MicroMagParameters RKSuitef90 ODE_Solvers LLODE_Debug util FortranCuda MagTenseMicroMagPyIO \
LandauLifshitzEquationSolver

# file names
LIBSRC_FILES = $(addsuffix .f90, ${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================
# names (without suffix), f90 sources
LIBSRC_WRAP_SOURCES = FortranToPythonIO

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90, ${LIBSRC_WRAP_SOURCES})

#=======================================================================
#                 Relevant suffixes
#=======================================================================
.SUFFIXES: .f90

#=======================================================================
#
#=======================================================================
.PHONY: all clean

all: ${PYTHON_MODN_ALL}

clean:
	${CLEAN_FILES}
	${CLEAN_FOLDER}
	${CLEAN_FOLDER_WIN}

.f90.o:
	${F90} ${F90FLAGS} -I $(OPT) -c $< -o $@

.c.o:
	${CC} ${CFLAGS} -c $< -o $@

${STATIC_LIB}: ${LIBSRC_OBJECTS} ${CUDA_OBJECTS}
	${LIBTOOL}

${PYTHON_MODN_ALL}: ${STATIC_LIB}
	f2py --build-dir . --fcompiler=${FCOMP} ${EXTRA_OPT} -c -m ${PYTHON_MODN} \
	${LIB_OPT} ${LIBSRC_WRAP_FILES} ${LIBSRC_OBJECTS} $(MKL) $(CUDA)
