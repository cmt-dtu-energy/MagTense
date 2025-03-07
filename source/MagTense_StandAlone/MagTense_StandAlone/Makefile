#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifx
FC = gfortran
FFLAGS =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = IO StandAloneIO MagTense_StandAlone
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	LIB_SUFFIX = .lib
else
	LIB_SUFFIX = .a
endif

NUM_INT_PATH = ../../NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../../TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = ../../DemagField/DemagField
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libTileDemagTensor${LIB_SUFFIX}
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: MagTense.x

clean:
	-rm -f *.x *.o *.mod

.f90.o:
	$(FC) $(FFLAGS) -I ${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} -c  $<

MagTense.x: ${LIBSRC_OBJECTS}
	$(FC) $(FFLAGS) StandAloneIO.o MagTense_StandAlone.o -Wl,--start-group \
		${DEMAG_FIELD_LIB} -Wl,--end-group -o MagTense.x
