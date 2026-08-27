#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifx
FFLAGS =
USE_CVODE =
CVODE_ROOT =
USE_MATLAB =
MATLAB_INCLUDE =
USE_CDFMM = 0
CDFMM_ROOT =
#=======================================================================
#       List all source files required for the project
#=======================================================================

LIBSRC_SOURCES = fmm3d_tree_mod MicroMagParameters fmm_nbor_tensor_mod util_micromag DemagAuxFunctions FortranCuda ComputeDifferentialOperatorsFromMesh UnstructuredMeshAnalysis

ifeq (${USE_CDFMM},1)
    CDFMM_FORTRAN_SOURCE = ${CDFMM_ROOT}/share/cdfmm/fortran/cdfmm_fortran.f90
    CDFMM_FORTRAN_OBJECT = cdfmm_fortran.o
    LIBSRC_SOURCES += DipFmmDemag
else
    CDFMM_FORTRAN_OBJECT =
endif

ifeq (${USE_MATLAB},1)
    OPT = -I${CONDA_PREFIX}/opt/compiler/include/intel64
	LIBSRC_SOURCES += MagTenseMicroMagIO
	MATLAB_OPT = -I${MATLAB_INCLUDE}
else
    OPT =
	LIBSRC_SOURCES += MagTenseMicroMagPyIO
	MATLAB_OPT =
endif


LIBSRC_SOURCES += LandauLifshitzEquationSolver
LIBSRC_OBJECTS = ${CDFMM_FORTRAN_OBJECT} $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	CONDA_PATH = $(subst \,/,${CONDA_PREFIX})
	MKL_ROOT = ${CONDA_PATH}/Library
	OPT = -I${CONDA_PATH}/opt/compiler/include/intel64
	LIBTOOL = LIB $^ /OUT:$@
	LIB_SUFFIX = .lib
else
	MKL_ROOT = ${CONDA_PREFIX}
	OPT =
	LIB_SUFFIX = .a

	ifeq (${UNAME}, Darwin)
		LIBTOOL = libtool -static -o $@ $^
	else
		LIBTOOL = ar src $@ $?
	endif
endif

STATIC_LIB = libMagTenseMicroMag${LIB_SUFFIX}
NUM_INT_PATH = ../NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../TileDemagTensor
DEMAG_FIELD_PATH = ../DemagField
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX}

AUXMT_PATH = ../AuxMT
AUXMT_LIB = ${AUXMT_PATH}/libAuxMT${LIB_SUFFIX}

#=======================================================================
#							Targets
#=======================================================================
# Clear the built-in suffix list before adding .f90. GNU make ships a .mod.o
# rule (Modula-2, m2c) and .mod comes before .f90 in the default list, so a
# leftover Fortran .mod file - which every interrupted build leaves behind -
# makes make build the object with m2c instead of the compiler.
.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: all clean

all: pre-build ${STATIC_LIB}
# Ignore .mod files as targets
%.o : %.mod

.f90.o:
	${FC} ${FFLAGS} -I${MKL_ROOT}/include -I${MKL_ROOT}/include/intel64/lp64 -I${AUXMT_PATH} -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} ${OPT} ${MATLAB_OPT} ${OPT} -c $< -o $@

ifeq (${USE_CDFMM},1)
cdfmm_fortran.o: ${CDFMM_FORTRAN_SOURCE}
	${FC} ${FFLAGS} -c $< -o $@

DipFmmDemag.o: cdfmm_fortran.o MicroMagParameters.o
LandauLifshitzEquationSolver.o: DipFmmDemag.o
endif

pre-build:
	cd ${NUM_INT_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"
	cd ${DEMAG_FIELD_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"

${STATIC_LIB}: ${DEMAG_FIELD_LIB} ${LIBSRC_OBJECTS}
	${LIBTOOL}

clean:
	-rm -f *${LIB_SUFFIX} *.o *.mod
