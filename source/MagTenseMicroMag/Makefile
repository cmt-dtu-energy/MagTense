#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifort

USE_MATLAB = 0
USE_CUDA = 1
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = MicroMagParameters LLODE_Debug FortranCuda \
MagTenseMicroMagPyIO LandauLifshitzEquationSolver

ifeq ($(USE_MATLAB),1)
	LIBSRC_SOURCES += MagTenseMicroMagIO
else
	LIBSRC_SOURCES += MagTenseMicroMagPyIO
endif

LIB_NAME = libMagTenseMicroMag

NUM_INT_LIB_NAME = libNumericalIntegration
NUM_INT_PATH = ../NumericalIntegration/NumericalIntegration

TILE_DEMAG_TENSOR_LIB_NAME = libTileDemagTensor
TILE_DEMAG_TENSOR_PATH = ../TileDemagTensor/TileDemagTensor

DEMAG_FIELD_LIB_NAME = libDemagField
DEMAG_FIELD_PATH = ../DemagField/DemagField

ifeq ($(OS),Windows_NT)
	FFLAGS = /nologo /real-size:64 /O3 /assume:nocc_omp /Qopenmp /fpp /fpe:0 \
	/fp:source /fp:precise /libs:static /DUSE_CUDA=${USE_CUDA} /DUSE_MATLAB=${USE_MATLAB}
	OPT = ${CONDA_PREFIX}/Library/include -I${MKL_ROOT}/include \
	-I${CONDA_PREFIX}/Library/include/intel64/ilp64
	F90_OPT = -I ${OPT}
	MKL_ROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
	MKL = -L${MKL_ROOT}/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64

	LIBTOOL = LIB $? /OUT:$@
	STATIC_LIB = ${LIB_NAME}.lib
	NUM_INT_LIB = ${NUM_INT_PATH}/${NUM_INT_LIB_NAME}.lib
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/${TILE_DEMAG_TENSOR_LIB_NAME}.lib
	DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/${DEMAG_FIELD_LIB_NAME}.lib
	CLEAN_FILES = -del ${STATIC_LIB} *.o *.mod

else
	FFLAGS = -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source \
	-fp-model precise -fpic -libs:static -DUSE_CUDA=${USE_CUDA} -DUSE_MATLAB=${USE_MATLAB}
	OPT = ${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/intel64/ilp64
	F90_OPT = -I ${OPT}
	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -lpthread -lm -ldl -liomp5 -lmkl_blas95_ilp64

	ifeq (${UNAME}, Darwin)
		LIBTOOL = libtool -static -o $@ $?
	else
		LIBTOOL = ar src $@ $?
	endif

	STATIC_LIB = ${LIB_NAME}.a
	NUM_INT_LIB = ${NUM_INT_PATH}/${NUM_INT_LIB_NAME}.a
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/${TILE_DEMAG_TENSOR_LIB_NAME}.a
	DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/${DEMAG_FIELD_LIB_NAME}.a
	CLEAN_FILES = -rm -f ${STATIC_LIB} *.o *.mod
endif

#=======================================================================
#						Relevant suffixes
#=======================================================================
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

.SUFFIXES: .f90

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: pre-build ${STATIC_LIB}

.f90.o:
	$(FC) $(FFLAGS) $(F90_OPT) -I ${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE)
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE)
	cd ${DEMAG_FIELD_PATH} && $(MAKE)

${STATIC_LIB}: ${LIBSRC_OBJECTS} ${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_LIB} ${DEMAG_FIELD_LIB}
	${LIBTOOL}

clean:
	${CLEAN_FILES}

