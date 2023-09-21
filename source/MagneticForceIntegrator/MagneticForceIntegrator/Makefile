#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran

MATLAB_INCLUDE = /apps/external/matlab/2019/b/extern/include/
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = MagForceIO MagForceMaxwellStressTensor MagneticForce

LIB_NAME = libMagneticForceIntegrator

NUM_INT_LIB_NAME = libNumericalIntegration
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration

TILE_DEMAG_TENSOR_LIB_NAME = libTileDemagTensor
TILE_DEMAG_TENSOR_PATH = ../../TileDemagTensor/TileDemagTensor

DEMAG_FIELD_LIB_NAME = libDemagField
DEMAG_FIELD_PATH = ../../DemagField/DemagField

OPT = -I ${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -I{MATLAB_INCLUDE}

ifeq ($(OS),Windows_NT)
	ifeq (${FC}, ifort)
		FFLAGS = /nologo /real-size:64 /O3 /assume:nocc_omp /Qopenmp /fpp \
		/fpe:0 /fp:source /fp:precise /libs:static /DUSE_CVODE=0
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 \
		-cpp -DUSE_CVODE=0
	endif

	LIBTOOL = LIB $? /OUT:$@
	STATIC_LIB = ${LIB_NAME}.lib
	NUM_INT_LIB = ${NUM_INT_LIB_NAME}.lib
	CLEAN_FILES = -del ${STATIC_LIB} *.o *.mod
else
	ifeq (${FC}, ifort)
		FFLAGS = -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source \
		-fp-model precise -fpic -libs:static -DUSE_CVODE=0
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 \
		-cpp -DUSE_CVODE=0
	endif

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	STATIC_LIB = ${LIB_NAME}.a
	NUM_INT_LIB = ${NUM_INT_LIB_NAME}.a
	CLEAN_FILES = -rm -f ${STATIC_LIB} *.o *.mod
endif

#=======================================================================
#						Relevant suffixes
#=======================================================================
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

.SUFFIXES: .for .f90 .F90 .f .o
#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: ${STATIC_LIB}

.f.o:
	$(FC) $(FFLAGS) -c $<
.f90.o:
	$(FC) $(FFLAGS) ${OPT} -c $<
.F90.o:
	$(FC) $(FFLAGS) -c $<
.for.o:
	$(FC) $(FFLAGS) -c $<

${STATIC_LIB}: ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	${CLEAN_FILES}

