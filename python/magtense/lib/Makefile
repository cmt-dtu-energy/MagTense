#=======================================================================
#                       compiler names and flags
#=======================================================================
USE_CUDA = 1
USE_CVODE = 0
USE_MICROMAG = 1

ifeq ($(OS),Windows_NT)
	F90 = ifort
else
	ifeq (${UNAME}, Darwin)
  		F90 = gfortran
	else
  		F90 = ifort
	endif
endif

ifeq (${F90}, gfortran)
	USE_CUDA = 0
	USE_MKL = 0
	USE_MICROMAG = 0
endif

CC = gcc
CFLAGS = -fPIC
CU = nvcc
PYTHON_MODN = magtensesource

VPATH = ../../../source/NumericalIntegration/NumericalIntegration:\
../../../source/TileDemagTensor/TileDemagTensor:\
../../../source/DemagField/DemagField:\
../../../source/MagTenseMicroMag

CU_PATH = ../../../source/MagTenseFortranCuda/cuda
CUDA_SOURCE = MagTenseCudaBlas
CUDA_WRAP_SOURCE = MagTenseCudaBlasICLWrapper

CUDA_FILE = $(addsuffix .cu, ${CUDA_SOURCE})
CUDA_WRAP_FILE = $(addsuffix .cxx, ${CUDA_WRAP_SOURCE})
CUDA_OBJECT = $(addsuffix .o, ${CUDA_SOURCE})
CUDA_WRAP = $(addsuffix .o, ${CUDA_WRAP_SOURCE})

ifeq ($(OS),Windows_NT)
	ifeq (${F90}, ifort)
		F90FLAGS = /nologo /O3 /assume:nocc_omp /Qopenmp /real-size:64 \
		/fp:precise /libs:static /threads /fpp /DUSE_CUDA=${USE_CUDA} \
		/DUSE_CVODE=${USE_CVODE} /DUSE_MICROMAG=${USE_MICROMAG} /DUSE_MATLAB=0
		FCOMP = intelvem
		EXTRA_OPT = -DNO_APPEND_FORTRAN -DUPPERCASE_FORTRAN --no-lower --opt='${F90FLAGS} /assume:nounderscore /names:uppercase'
		OPT = ${CONDA_PREFIX}/Library/include -I${CONDA_PREFIX}/Library/include/intel64/ilp64
		F90_OPT = -I ${OPT}
		LIB_OPT = -I:${OPT} -L. -llibsrc
	else ifeq (${F90}, gfortran)
		F90FLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp -DUSE_MICROMAG=${USE_MICROMAG}
		FCOMP = gnu95 --compiler=mingw32
		EXTRA_OPT = --opt='${F90FLAGS}'
		F90_OPT =
		LIB_OPT = -L. -llibsrc
	endif

	CPP = icl
	CPPFLAGS = /Fo:${CUDA_WRAP}
	CUFLAGS = -I ${OPT}
	CUDA = -L${CONDA_PREFIX}/Lib/x64 -lcublas -lcudart -lcuda -lcusparse
	MKL = -L${CONDA_PREFIX}/Library/lib -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
	MKL_OMP = -L"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\compiler\lib\intel64_win" -llibiomp5md
		
	LIBTOOL = LIB $? /OUT:$@
	STATIC_LIB = libsrc.lib
	CLEAN_FILES = -del libsrc.lib *.o *.mod *.pyd
	CLEAN_FOLDER = -rmdir /s /q Release 2>NUL || VER>NUL
	CLEAN_FOLDER_WIN = -for /f %%i in ('dir /a:d /b src.win-*') do rd /s /q %%i
	PYTHON_MODN_ALL = _${PYTHON_MODN}.dll

else
	export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib
	
	ifeq (${F90}, ifort)
		F90FLAGS = -fPIC -nologo -O3 -assume nocc_omp -qopenmp -real-size 64 \
		-fp-model precise -threads -fpp -DUSE_CUDA=${USE_CUDA} -DUSE_CVODE=${USE_CVODE} \
		-DUSE_MICROMAG=${USE_MICROMAG} -DUSE_MATLAB=0
		FCOMP = intelem
		OPT = ${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/intel64/ilp64
		F90_OPT = -I ${OPT}
		LIB_OPT = -I${OPT} -L. -lsrc
	else ifeq (${F90}, gfortran)
		F90FLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp -DUSE_MICROMAG=${USE_MICROMAG}
		FCOMP = gnu95
		F90_OPT =
		LIB_OPT = -L. -lsrc
	endif

	EXTRA_OPT = --opt='${F90FLAGS}'
	CPP = icx
	CPPFLAGS = ${CFLAGS} -o ${CUDA_WRAP}
	CUFLAGS = -shared -Xcompiler -fPIC -I ${OPT}
	CUDA = -L${CONDA_PREFIX}/lib -lcublas -lcudart -lcuda -lcusparse
	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -lpthread -lm -ldl -liomp5 -lmkl_blas95_lp64

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

ifeq ($(USE_MKL),0)
	MKL =
endif

ifeq ($(USE_CUDA),0)
	CUDA =
	CUDA_WRAP =
	CUDA_OBJECT =
endif

#=======================================================================
#       List all source files required for the project
#=======================================================================
# names (without suffix), f90 sources
LIBSRC_SOURCES = IOgeneral IntegrationDataTypes quadpack SpecialFunctions TileTensorHelperFunctions \
TileRectangularPrismTensor TileCircPieceTensor TileCylPieceTensor TilePlanarCoilTensor \
TileTriangle TileNComponents DemagFieldGetSolution spline MagParameters IterateMagnetSolution

ifeq ($(USE_MICROMAG),1)
	LIBSRC_SOURCES += MicroMagParameters RKSuitef90 ODE_Solvers LLODE_Debug util FortranCuda \
	MagTenseMicroMagPyIO LandauLifshitzEquationSolver
endif

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
	${F90} ${F90FLAGS} $(F90_OPT) -c $< -o $@

.c.o:
	${CC} ${CFLAGS} -c $< -o $@

${CUDA_OBJECT}:
	${CU} ${CUFLAGS} -c ${CU_PATH}/${CUDA_FILE} -o $@

${CUDA_WRAP}:
	${CPP} ${CPPFLAGS} -c ${CU_PATH}/${CUDA_WRAP_FILE}

${STATIC_LIB}: ${LIBSRC_OBJECTS} ${CUDA_OBJECT} ${CUDA_WRAP}
	${LIBTOOL}

${PYTHON_MODN_ALL}: ${STATIC_LIB}
	f2py --build-dir . --fcompiler=${FCOMP} ${EXTRA_OPT} -c -m ${PYTHON_MODN} \
	${LIB_OPT} ${LIBSRC_WRAP_FILES} ${LIBSRC_OBJECTS} $(MKL) $(CUDA)
