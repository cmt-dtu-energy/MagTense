#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifort
USE_MATLAB = 0
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

LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
OPT = ${CONDA_PREFIX}/Library/include -I${MKL_ROOT}/include -I${CONDA_PREFIX}/Library/include/intel64/ilp64
F90_OPT = -I ${OPT}
MKL_ROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
MKL = -L${MKL_ROOT}/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64

LIBTOOL = LIB $? /OUT:$@
LIB_SUFFIX = .lib
CLEAN_FILES = -del *.lib *.o *.mod

else
OPT = ${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/intel64/ilp64
F90_OPT = -I ${OPT}
MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -lpthread -lm -ldl -liomp5 -lmkl_blas95_ilp64

LIB_SUFFIX = .a
CLEAN_FILES = -rm -f *.a *.o *.mod

ifeq (${UNAME}, Darwin)
	LIBTOOL = libtool -static -o $@ $?
else
	LIBTOOL = ar src $@ $?
endif
endif

STATIC_LIB = libMagTenseMicroMag${LIB_SUFFIX}
NUM_INT_PATH = ../NumericalIntegration/NumericalIntegration
NUM_INT_LIB = ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX}

TILE_DEMAG_TENSOR_PATH = ../TileDemagTensor/TileDemagTensor
TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX}

DEMAG_FIELD_PATH = ../DemagField/DemagField
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX}
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: pre-build ${STATIC_LIB}

.f90.o:
	$(FC) $(FFLAGS) $(F90_OPT) -I ${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd ${DEMAG_FIELD_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'

${STATIC_LIB}: ${LIBSRC_OBJECTS} ${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_LIB} ${DEMAG_FIELD_LIB}
	${LIBTOOL}

clean:
	${CLEAN_FILES}

