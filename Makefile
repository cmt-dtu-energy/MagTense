#=======================================================================
#                       compiler names and flags
#=======================================================================
USE_CUDA = 1
USE_CVODE = 0
USE_MICROMAG = 1
USE_MATLAB = 0
MATLAB_INCLUDE =
# /usr/local/MATLAB/<version>/extern/include (Linux)
# "C:\Program Files\MATLAB\<version>\extern\include" (Win)
CPP = icx
FC = ifx
MKFILE_PATH := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

ifeq (${UNAME}, Darwin)
	FC = gfortran
	USE_CUDA = 0
	USE_MKL = 0
	USE_MICROMAG = 0
	USE_MATLAB = 0
endif

ifeq (${FC}, ifx)
	ifeq ($(OS),Windows_NT)
		FFLAGS = /O3 /fpp /real-size:64 /Qopenmp /assume:nocc_omp /fpe:0 \
			/fp:source /nologo /DUSE_CVODE=${USE_CVODE} /DUSE_MATLAB=${USE_MATLAB} \
			/DUSE_CUDA=${USE_CUDA} /DUSE_MICROMAG=${USE_MICROMAG}
	else
		FFLAGS = -O3 -fpp -real-size 64 -qopenmp -assume nocc_omp -fpe0 \
			-fp-model=source -fpic -nologo -DUSE_CVODE=${USE_CVODE} \
			-DUSE_MATLAB=${USE_MATLAB} -DUSE_CUDA=${USE_CUDA} \
			-DUSE_MICROMAG=${USE_MICROMAG}
	endif
else ifeq (${FC}, gfortran)
	FFLAGS = -O3 -fdefault-real-8 -fopenmp -ffree-line-length-512 -cpp -fPIC \
		-DUSE_MICROMAG=0 -DUSE_CVODE=${USE_CVODE}
endif

PYTHON_MODN = magtensesource
PYTHON_LIBPATH = python/src/magtense/lib

NUM_INT_PATH = source/NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = source/TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = source/DemagField/DemagField
MICROMAG_PATH = source/MagTenseMicroMag
FORTRAN_CUDA_PATH = source/MagTenseFortranCuda/cuda
STANDALONE_PATH = source/MagTense_StandAlone/MagTense_StandAlone
FORCEINTEGRATOR_PATH = source/MagneticForceIntegrator/MagneticForceIntegrator

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
	PY_MOD_SUFFIX = .pyd

	ifeq (${FC}, ifx)
		EXTRA_FFLAGS = "${FFLAGS} /assume:underscore /names:lowercase"
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
	PY_MOD_SUFFIX = .so

	ifeq (${FC}, ifx)
		EXTRA_FFLAGS = "${FFLAGS}"
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

ifeq ($(USE_MATLAB),0)
	FORCEINTEGRATOR =
else
	FORCEINTEGRATOR = forceintegrator
	CP_LIB += && cp ${FORCEINTEGRATOR_PATH}/libMagneticForceIntegrator${LIB_SUFFIX} .
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

PYTHON_MODN_ALL = _${PYTHON_MODN}${PY_MOD_SUFFIX}

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR}

standalone: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} standalone

python: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${PYTHON_MODN_ALL}

ps: magnetostatic ${MICROMAG}

cmdx64: ${PYTHON_MODN_ALL}

clean:
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	cd $(STANDALONE_PATH) && ${MAKE} clean
	cd $(FORCEINTEGRATOR_PATH) && ${MAKE} clean
	rm -f *${LIB_SUFFIX} *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -f -r ${PYTHON_LIBPATH}/build

magnetostatic:
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd ${DEMAG_FIELD_PATH}  && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)

micromagnetism:
	cd ${MICROMAG_PATH} && $(MAKE) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)

cuda:
	cd ${FORTRAN_CUDA_PATH} && $(MAKE) CPP=$(CPP)

forceintegrator:
	cd $(FORCEINTEGRATOR_PATH) && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' MATLAB_INCLUDE=$(MATLAB_INCLUDE)

standalone:
	cd $(STANDALONE_PATH) && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}'
	mkdir build
	cp $(STANDALONE_PATH)/MagTense.x build/MagTense.x

${PYTHON_MODN_ALL}:
	${CP_LIB}
	FC=${FC} FFLAGS=${EXTRA_FFLAGS} LDFLAGS=${LDFLAGS} \
		python -m numpy.f2py -c -m ${PYTHON_MODN} \
		--build-dir ${PYTHON_LIBPATH}/build -I${OPT} -I${INCLUDE_OBJ} ${CVODE_OPT} \
		-L${MKFILE_PATH} ${LIB_OPT} python/FortranToPythonIO.f90 ${MKL} ${CUDA} ${CVODE}
	cp *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/
