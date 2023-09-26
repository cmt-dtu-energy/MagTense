#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran

#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = IO StandAloneIO MagTense_StandAlone

NUM_INT_LIB_NAME = libNumericalIntegration
NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration

TILE_DEMAG_TENSOR_LIB_NAME = libTileDemagTensor
TILE_DEMAG_TENSOR_PATH = ../../TileDemagTensor/TileDemagTensor

DEMAG_FIELD_LIB_NAME = libDemagField
DEMAG_FIELD_PATH = ../../DemagField/DemagField

ifeq ($(OS),Windows_NT)
	ifeq (${FC}, ifort)
		FFLAGS = /nologo /real-size:64 /O3 /assume:nocc_omp /Qopenmp /fpp \
		/fpe:0 /fp:source /fp:precise
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp
	endif

	NUM_INT_LIB = ${NUM_INT_LIB_NAME}.lib
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_LIB_NAME}.lib
	DEMAG_FIELD_LIB = ${DEMAG_FIELD_LIB_NAME}.lib
	CLEAN_FILES = -del *.x *.o *.mod
else
	ifeq (${FC}, ifort)
		FFLAGS = -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source \
		-fp-model precise -fpic
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp
	endif

	NUM_INT_LIB = ${NUM_INT_LIB_NAME}.a
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_LIB_NAME}.a
	DEMAG_FIELD_LIB = ${DEMAG_FIELD_LIB_NAME}.a
	CLEAN_FILES = -rm -f *.x *.o *.mod
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

all: MagTense.x

clean:
	${CLEAN_FILES}

.f90.o:
	$(FC) $(FFLAGS) -I ${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -c  $<

MagTense.x: ${LIBSRC_OBJECTS}
	$(FC) $(FFLAGS) StandAloneIO.o MagTense_StandAlone.o -Wl,--start-group \
	${NUM_INT_PATH}/${NUM_INT_LIB} \
	${TILE_DEMAG_TENSOR_PATH}/${TILE_DEMAG_TENSOR_LIB} \
	${DEMAG_FIELD_PATH}/${DEMAG_FIELD_LIB} \
	-Wl,--end-group -o MagTense.x
