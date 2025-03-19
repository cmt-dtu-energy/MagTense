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
LIBSRC_SOURCES = MicroMagParameters util_micromag DemagAuxFunctions FortranCuda

ifeq (${USE_MATLAB},1)
	LIBSRC_SOURCES += MagTenseMicroMagIO
	MATLAB_OPT = -I"${MATLAB_INCLUDE}"
else
	LIBSRC_SOURCES += MagTenseMicroMagPyIO
	MATLAB_OPT =
endif

LIBSRC_SOURCES += LandauLifshitzEquationSolver
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	CONDA_PATH = $(subst \,/,${CONDA_PREFIX})
	MKL_ROOT = ${CONDA_PATH}/Library
	OPT = -I${CONDA_PATH}/opt/compiler/include/intel64
	LIBTOOL = LIB $? /NOLOGO /OUT:$@
	LIB_SUFFIX = .lib
else
	MKL_ROOT = ${CONDA_PREFIX}
	OPT =
	LIB_SUFFIX = .a

	ifeq (${UNAME}, Darwin)
		LIBTOOL = libtool -static -o $@ $?
	else
		LIBTOOL = ar src $@ $?
	endif
endif

STATIC_LIB = libMagTenseMicroMag${LIB_SUFFIX}
NUM_INT_PATH = ../NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = ../DemagField/DemagField
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX}

#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: pre-build ${STATIC_LIB}

.f90.o:
	${FC} ${FFLAGS} -I${MKL_ROOT}/include -I${MKL_ROOT}/include/intel64/lp64 -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} ${OPT} ${MATLAB_OPT} -c $< -o $@

pre-build:
	cd ${NUM_INT_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"
	cd ${DEMAG_FIELD_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"

${STATIC_LIB}: ${DEMAG_FIELD_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	-rm -f *${LIB_SUFFIX} *.o *.mod
