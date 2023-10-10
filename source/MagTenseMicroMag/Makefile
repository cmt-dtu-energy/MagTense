#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifort
USE_MATLAB = 0
MATLAB_INCLUDE =
MKL_ROOT =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = MicroMagParameters LLODE_Debug FortranCuda

ifeq ($(USE_MATLAB),1)
	LIBSRC_SOURCES += MagTenseMicroMagIO
else
	LIBSRC_SOURCES += MagTenseMicroMagPyIO
	MATLAB_INCLUDE =
endif

LIBSRC_SOURCES += LandauLifshitzEquationSolver
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
LIBTOOL = LIB $? /OUT:$@
LIB_SUFFIX = .lib
CLEAN_FILES = -del *.lib *.o *.mod

else
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
	$(FC) $(FFLAGS) -I${MKL_ROOT}/include -I${MKL_ROOT}/include/intel64/ilp64 -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -I${MATLAB_INCLUDE} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd ${DEMAG_FIELD_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'

${STATIC_LIB}: ${LIBSRC_OBJECTS} ${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_LIB} ${DEMAG_FIELD_LIB}
	${LIBTOOL}

clean:
	${CLEAN_FILES}

