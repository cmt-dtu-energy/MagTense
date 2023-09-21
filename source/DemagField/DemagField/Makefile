#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran

#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = DemagFieldGetSolution MagParameters spline IterateMagnetSolution 

LIB_NAME = libDemagField

NUM_INT_LIB_NAME = libNumericalIntegration
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration

TILE_DEMAG_TENSOR_LIB_NAME = libTileDemagTensor
TILE_DEMAG_TENSOR_PATH = ../../TileDemagTensor/TileDemagTensor

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
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_LIB_NAME}.lib
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
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_LIB_NAME}.a
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

all: pre-build ${STATIC_LIB}

.f.o:
	$(FC) $(FFLAGS) -c $< -o $@
.f90.o:
	$(FC) $(FFLAGS) -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -c $< -o $@
.F90.o:
	$(FC) $(FFLAGS) -c $< -o $@
.for.o:
	$(FC) $(FFLAGS) -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE)
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE)

${STATIC_LIB}: ${NUM_INT_PATH}/${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_PATH}/${TILE_DEMAG_TENSOR_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	${CLEAN_FILES}
