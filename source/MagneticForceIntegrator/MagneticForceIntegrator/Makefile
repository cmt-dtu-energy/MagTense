#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran
FFLAGS =

MATLAB_INCLUDE =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = MagForceIO MagForceMaxwellStressTensor MagneticForce
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
LIB_SUFFIX = .lib
LIBTOOL = LIB $? /OUT:$@
CLEAN_FILES = -del ${STATIC_LIB} *.o *.mod
else
ifeq (${UNAME}, Darwin)
	LIBTOOL = libtool -static -o $@ $?
else
	LIBTOOL = ar src $@ $?
endif
LIB_SUFFIX = .a
CLEAN_FILES = -rm -f ${STATIC_LIB} *.o *.mod
endif

STATIC_LIB = libMagneticForceIntegrator${LIB_SUFFIX}
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration
NUM_INT_LIB = ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX}

TILE_DEMAG_TENSOR_PATH = ../../TileDemagTensor/TileDemagTensor
TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX}

DEMAG_FIELD_PATH = ../../DemagField/DemagField
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX}
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: ${STATIC_LIB}

.f90.o:
	$(FC) $(FFLAGS) -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -I${MATLAB_INCLUDE} -c $<

${STATIC_LIB}: ${LIBSRC_OBJECTS} ${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_LIB} ${DEMAG_FIELD_LIB}
	${LIBTOOL}

clean:
	${CLEAN_FILES}
