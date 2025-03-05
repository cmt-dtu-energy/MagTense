#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifx
FC = gfortran
FFLAGS =
USE_CVODE =
CVODE_ROOT =
USE_MATLAB =
MATLAB_INCLUDE =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = DemagFieldGetSolution MagParameters spline IterateMagnetSolution

ifeq ($(USE_MATLAB),1)
	LIBSRC_SOURCES += MagTileIO
	MATLAB_OPT = -I"${MATLAB_INCLUDE}"
else
	MATLAB_OPT =
endif

LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	OPT = -I${CONDA_PREFIX}/opt/compiler/include/intel64
	LIBTOOL = LIB $? /NOLOGO /OUT:$@
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

STATIC_LIB = libDemagField${LIB_SUFFIX}
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration
NUM_INT_LIB = ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX}

TILE_DEMAG_TENSOR_PATH = ../../TileDemagTensor/TileDemagTensor
TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX}
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: pre-build ${STATIC_LIB}

.f90.o:
	$(FC) $(FFLAGS) -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} ${OPT} ${MATLAB_OPT} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT="$(CVODE_ROOT)" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="$(MATLAB_INCLUDE)"
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT="$(CVODE_ROOT)" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="$(MATLAB_INCLUDE)"

${STATIC_LIB}: ${TILE_DEMAG_TENSOR_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	-rm -f *${LIB_SUFFIX} *.o *.mod
