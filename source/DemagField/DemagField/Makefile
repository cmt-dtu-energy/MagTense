#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran
FFLAGS =
USE_MATLAB = 0
MATLAB_INCLUDE =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = DemagFieldGetSolution MagParameters spline IterateMagnetSolution

ifeq ($(USE_MATLAB),1)
	LIBSRC_SOURCES += MagTileIO
endif

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
	$(FC) $(FFLAGS) -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${MATLAB_INCLUDE} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'

${STATIC_LIB}: ${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	${CLEAN_FILES}
