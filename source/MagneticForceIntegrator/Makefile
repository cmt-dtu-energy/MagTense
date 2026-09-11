#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifx
FFLAGS =
MATLAB_INCLUDE =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = MagForceIO MagForceMaxwellStressTensor MagneticForce
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	OPT = -I${CONDA_PREFIX}/opt/compiler/include/intel64
	LIB_SUFFIX = .lib
	LIBTOOL = LIB $^ /OUT:$@
else
	OPT =
	LIB_SUFFIX = .a

	ifeq (${UNAME}, Darwin)
		LIBTOOL = libtool -static -o $@ $^
	else
		LIBTOOL = ar src $@ $?
	endif
endif

STATIC_LIB = libMagneticForceIntegrator${LIB_SUFFIX}
AUXMT_PATH = ../AuxMT
NUM_INT_PATH = ../NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../TileDemagTensor
DEMAG_FIELD_PATH = ../DemagField
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX}
#=======================================================================
#							Targets
#=======================================================================
# Clear the built-in suffix list before adding .f90. GNU make ships a .mod.o
# rule (Modula-2, m2c) and .mod comes before .f90 in the default list, so a
# leftover Fortran .mod file - which every interrupted build leaves behind -
# makes make build the object with m2c instead of the compiler.
.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: all clean

all: ${STATIC_LIB}

.f90.o:
	$(FC) $(FFLAGS) -I${AUXMT_PATH} -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -I${MATLAB_INCLUDE} ${OPT} -c $< -o $@

${STATIC_LIB}: ${DEMAG_FIELD_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	-rm -f ${STATIC_LIB} *.o *.mod
