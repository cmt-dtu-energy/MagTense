#==========================
# Compiler names and flags
#==========================
USE_CUDA = 1
USE_CVODE = 0
USE_MICROMAG = 1
CPP = icx
FC = ifx
MKFILE_PATH := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

ifeq (${UNAME}, Darwin)
	FC = gfortran
	USE_CUDA = 0
	USE_MKL = 0
	USE_MICROMAG = 0
endif

ifeq (${FC}, ifx)
	ifeq ($(OS),Windows_NT)
		FFLAGS = /O3 /fpp /real-size:64 /Qopenmp /assume:nocc_omp /fpe:0 \
			/fp:source /nologo /DUSE_CVODE=${USE_CVODE} /DUSE_MATLAB=0 \
			/DUSE_CUDA=${USE_CUDA} /DUSE_MICROMAG=${USE_MICROMAG}
	else
		FFLAGS = -O3 -fpp -real-size 64 -qopenmp -assume nocc_omp -fpe0 \
			-fp-model=source -fpic -nologo -DUSE_CVODE=${USE_CVODE} \
			-DUSE_MATLAB=0 -DUSE_CUDA=${USE_CUDA} -DUSE_MICROMAG=${USE_MICROMAG}
	endif
else ifeq (${FC}, gfortran)
	FFLAGS = -O3 -fdefault-real-8 -fopenmp -ffree-line-length-512 -cpp -fPIC \
		-DUSE_MICROMAG=0 -DUSE_CVODE=${USE_CVODE}
endif

PYTHON_MODN = magtensesource

NUM_INT_PATH = ../../../../source/NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../../../../source/TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = ../../../../source/DemagField/DemagField
MICROMAG_PATH = ../../../../source/MagTenseMicroMag
FORTRAN_CUDA_PATH = ../../../../source/MagTenseFortranCuda/cuda

VPATH = ${NUM_INT_PATH}:${TILE_DEMAG_TENSOR_PATH}:${DEMAG_FIELD_PATH}:\
${MICROMAG_PATH}:${FORTRAN_CUDA_PATH}

ifeq ($(OS),Windows_NT)
	CONDA_PATH = $(subst \,/,${CONDA_PREFIX})
	CUDA = -L${CONDA_PATH}/Library/lib -lcublas -lcudart -lcuda -lcusparse
	MKL = -L${CONDA_PATH}/Library/lib -lmkl_intel_lp64_dll -lmkl_intel_thread_dll \
		-lmkl_core_dll -lmkl_blas95_lp64 -llibiomp5md
	CVODE_ROOT = "C:\Program Files (x86)\sundials-4.1.0\instdir"
	LDFLAGS = '/DEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcmt.lib'
	LIB_SUFFIX = .lib
	CLEAN_FILES = -rm *.lib *.o *.mod *.pyd
	CLEAN_FOLDER = -rm -r build
	PYTHON_MODN_ALL = _${PYTHON_MODN}.dll

	ifeq (${FC}, ifx)
		EXTRA_FFLAGS = '${FFLAGS} /assume:underscore /names:lowercase'
		OPT = ${CONDA_PATH}/Library/include \
			-I${CONDA_PATH}/Library/include/intel64/lp64 \
			-I${CONDA_PATH}/opt/compiler/include/intel64
		LIB_OPT = -llibDemagField
	endif
else
	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -liomp5 -lmkl_blas95_lp64 -lpthread -lm -ldl
	CUDA = -L${CONDA_PREFIX}/lib -lcublas -lcudart -lcuda -lcusparse
	CVODE_ROOT = "/usr/local/sundials-4.1.0/instdir"
	LDFLAGS =
	LIB_SUFFIX = .a
	CLEAN_FILES = -rm *.a *.o *.mod *.so 
	CLEAN_FOLDER = -rm -r build
	PYTHON_MODN_ALL = _${PYTHON_MODN}.so

	ifeq (${FC}, ifx)
		EXTRA_FFLAGS = '${FFLAGS}'
		OPT = ${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/intel64/lp64
		LIB_OPT = -lNumericalIntegration -lTileDemagTensor -lDemagField
	endif
endif

ifeq (${FC}, gfortran)
	LIB_OPT =
endif

ifeq ($(USE_MKL),0)
	MKL =
endif

ifeq (${USE_CVODE},0)
	CVODE =
	CVODE_OPT =
else
	CVODE = -L${CVODE_ROOT}/lib -lsundials_nvecserial -lsundials_sunmatrixdense \
		-lsundials_sunlinsoldense -lsundials_fnvecserial_mod -lsundials_cvode \
		-lsundials_fsunnonlinsolfixedpoint_mod
	CVODE_OPT = -I${CVODE_ROOT}/fortran
endif

ifeq ($(USE_MICROMAG),0)
	MICROMAG =
	ifeq ($(OS),Windows_NT)
		CP_LIB = cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
	else
		CP_LIB = cp ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX} .
		CP_LIB += && cp ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX} .
		CP_LIB += && cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
	endif
else
	MICROMAG = micromagnetism
	ifeq ($(OS),Windows_NT)
		LIB_OPT = -llibMagTenseMicroMag
		CP_LIB = cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
	else
		LIB_OPT += -lMagTenseMicroMag
		CP_LIB = cp ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX} .
		CP_LIB += && cp ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX} .
		CP_LIB += && cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
		CP_LIB += && cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
	endif
endif

ifeq ($(USE_CUDA),0)
	COMPILE_CUDA =
	CUDA =
else
	COMPILE_CUDA = cuda
	CP_LIB += && cp ${FORTRAN_CUDA_PATH}/libCuda${LIB_SUFFIX} .
	ifeq ($(OS),Windows_NT)
		LIB_OPT += -llibCuda
	else
		LIB_OPT += -lCuda
	endif
endif

INCLUDE_OBJ = ${MKFILE_PATH}/${NUM_INT_PATH} \
	-I${MKFILE_PATH}/${TILE_DEMAG_TENSOR_PATH} \
	-I${MKFILE_PATH}/${DEMAG_FIELD_PATH} \
	-I${MKFILE_PATH}/${MICROMAG_PATH} \
	-I${MKFILE_PATH}/${FORTRAN_CUDA_PATH}

LIBSRC_SOURCES = ${NUM_INT_PATH}/IOgeneral \
	${NUM_INT_PATH}/IntegrationDataTypes \
	${NUM_INT_PATH}/quadpack \
	${NUM_INT_PATH}/SpecialFunctions \
	${TILE_DEMAG_TENSOR_PATH}/TileTensorHelperFunctions \
	${TILE_DEMAG_TENSOR_PATH}/TileRectangularPrismTensor \
	${TILE_DEMAG_TENSOR_PATH}/TileCircPieceTensor \
	${TILE_DEMAG_TENSOR_PATH}/TileCylPieceTensor \
	${TILE_DEMAG_TENSOR_PATH}/TilePlanarCoilTensor \
	${TILE_DEMAG_TENSOR_PATH}/TileTriangle \
	${TILE_DEMAG_TENSOR_PATH}/TileNComponents \
	${DEMAG_FIELD_PATH}/DemagFieldGetSolution \
	${DEMAG_FIELD_PATH}/spline \
	${DEMAG_FIELD_PATH}/MagParameters \
	${DEMAG_FIELD_PATH}/IterateMagnetSolution

ifeq ($(USE_MICROMAG),1)
	LIBSRC_SOURCES += ${MICROMAG_PATH}/MicroMagParameters \
		${NUM_INT_PATH}/RKSuitef90 \
		${NUM_INT_PATH}/ODE_Solvers \
		${MICROMAG_PATH}/LLODE_Debug \
		${NUM_INT_PATH}/util \
		${MICROMAG_PATH}/FortranCuda \
		${MICROMAG_PATH}/MagTenseMicroMagPyIO \
		${MICROMAG_PATH}/LandauLifshitzEquationSolver
endif

LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all ps cmdx64 clean clean_build

all: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${PYTHON_MODN_ALL}

ps: magnetostatic ${MICROMAG}

cmdx64: ${PYTHON_MODN_ALL}

clean:
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	${CLEAN_FILES}
	${CLEAN_FOLDER}

clean_build:
	-rm *.pyd
	${CLEAN_FOLDER}

magnetostatic:
	cd ${NUM_INT_PATH} && ${MAKE} FC=${FC} FFLAGS='${FFLAGS}' CVODE_OPT=${CVODE_OPT} USE_MATLAB=0
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} FC=${FC} FFLAGS='${FFLAGS}'
	cd ${DEMAG_FIELD_PATH} && ${MAKE} FC=${FC} FFLAGS='${FFLAGS}' USE_MATLAB=0

micromagnetism:
	cd ${MICROMAG_PATH} && ${MAKE} FC=${FC} FFLAGS='${FFLAGS}' USE_CVODE=${USE_CVODE} CVODE_ROOT=${CVODE_ROOT} USE_MATLAB=0 MATLAB_INCLUDE=""

cuda:
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} CPP=${CPP}

${PYTHON_MODN_ALL}:
	${CP_LIB}
	FC=${FC} FFLAGS=${EXTRA_FFLAGS} LDFLAGS=${LDFLAGS} python -m numpy.f2py -c -m ${PYTHON_MODN} \
		--build-dir ./build -I${OPT} -I${INCLUDE_OBJ} ${CVODE_OPT} \
		-L${MKFILE_PATH} ${LIB_OPT} FortranToPythonIO.f90 ${MKL} ${CUDA} ${CVODE}
