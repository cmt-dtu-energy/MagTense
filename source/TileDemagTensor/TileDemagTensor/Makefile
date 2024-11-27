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
LIBSRC_SOURCES = TileTensorHelperFunctions TileCylPieceTensor \
TileRectangularPrismTensor TilePlanarCoilTensor TileTriangle \
TileCircPieceTensor TileNComponents
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	LIBTOOL = LIB $? /OUT:$@
	LIB_SUFFIX = .lib
	CLEAN_FILES = -del *.lib *.o *.mod
else
	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	LIB_SUFFIX = .a
	CLEAN_FILES = -rm -f *.a *.o *.mod
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
	$(FC) $(FFLAGS) -I ${NUM_INT_PATH} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT="$(CVODE_ROOT)" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="$(MATLAB_INCLUDE)"

${STATIC_LIB}: ${NUM_INT_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	${CLEAN_FILES}
	