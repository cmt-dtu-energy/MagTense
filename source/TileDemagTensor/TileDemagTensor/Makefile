#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran

#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = TileTensorHelperFunctions TileCylPieceTensor \
TileRectangularPrismTensor TilePlanarCoilTensor TileTriangle \
TileCircPieceTensor TileNComponents

LIB_NAME = libTileDemagTensor

NUM_INT_LIB_NAME = libNumericalIntegration
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration

ifeq ($(OS),Windows_NT)
	ifeq (${FC}, ifort)
		FFLAGS = /nologo /real-size:64 /O3 /assume:nocc_omp /Qopenmp /fpp \
		/fpe:0 /fp:source /fp:precise /libs:static
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp
	endif

	LIBTOOL = LIB $? /OUT:$@
	STATIC_LIB = ${LIB_NAME}.lib
	NUM_INT_LIB = ${NUM_INT_PATH}/${NUM_INT_LIB_NAME}.lib
	CLEAN_FILES = -del ${STATIC_LIB} *.o *.mod
else
	ifeq (${FC}, ifort)
		FFLAGS = -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source \
		-fp-model precise -fpic -libs:static
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp
	endif

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	STATIC_LIB = ${LIB_NAME}.a
	NUM_INT_LIB = ${NUM_INT_PATH}/${NUM_INT_LIB_NAME}.a
	CLEAN_FILES = -rm -f ${STATIC_LIB} *.o *.mod
endif

#=======================================================================
#						Relevant suffixes
#=======================================================================
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

.SUFFIXES: .f90

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: pre-build ${STATIC_LIB}

.f90.o:
	$(FC) $(FFLAGS) -I ${NUM_INT_PATH} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && $(MAKE)

${STATIC_LIB}: ${NUM_INT_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	${CLEAN_FILES}
	