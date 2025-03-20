#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifx
FFLAGS =
USE_CVODE =
CVODE_ROOT =
USE_MATLAB =
MATLAB_INCLUDE =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = TileTensorHelperFunctions TileCylPieceTensor \
	TileRectangularPrismTensor TilePlanarCoilTensor TileTriangle \
	TileCircPieceTensor TileNComponents
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	OPT = -I${CONDA_PREFIX}/opt/compiler/include/intel64
	LIBTOOL = LIB $? /OUT:$@
	LIB_SUFFIX = .lib
else
	OPT =
	LIB_SUFFIX = .a

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
endif

STATIC_LIB = libTileDemagTensor${LIB_SUFFIX}
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration
NUM_INT_LIB = ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX}
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: pre-build ${STATIC_LIB}

.f90.o:
	${FC} ${FFLAGS} -I ${NUM_INT_PATH} ${OPT} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && ${MAKE} FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="${MATLAB_INCLUDE}"

${STATIC_LIB}: ${NUM_INT_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	-rm -f *${LIB_SUFFIX} *.o *.mod
	