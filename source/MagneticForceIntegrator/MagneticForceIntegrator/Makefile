#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifx
FC = gfortran
FFLAGS =
MATLAB_INCLUDE =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = MagForceIO MagForceMaxwellStressTensor MagneticForce

ifeq ($(OS),Windows_NT)
	OPT = -I${CONDA_PREFIX}/opt/compiler/include/intel64
	LIB_SUFFIX = .lib
	LIBTOOL = LIB $? /OUT:$@
else
	OPT =
	LIB_SUFFIX = .a

	ifeq (${UNAME}, Darwin)
		LIBTOOL = libtool -static -o $@ $?
	else
		LIBTOOL = ar src $@ $?
	endif
endif

LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})
STATIC_LIB = libMagneticForceIntegrator${LIB_SUFFIX}
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../../TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = ../../DemagField/DemagField
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX}
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: ${STATIC_LIB}

.f90.o:
	$(FC) $(FFLAGS) -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -I${MATLAB_INCLUDE} ${OPT} -c $< -o $@

${STATIC_LIB}: ${DEMAG_FIELD_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	-rm -f ${STATIC_LIB} *.o *.mod
