#=======================================================================
#                       compiler names and flags
#=======================================================================
USE_CUDA = 1
USE_CVODE = 0
USE_MICROMAG = 1
USE_MATLAB = 0
MATLAB_INCLUDE =
USE_FMM3D ?= 0

# /usr/local/MATLAB/<version>/extern/include (Linux)
# "C:\Program Files\MATLAB\<version>\extern\include" (Win)
CPP = icx
FC = ifx
MKFILE_PATH := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
CVODE_ROOT = ${MKFILE_PATH}/cvode

#=======================================================================
#                    FMM3D integration (upstream Makefile)
#=======================================================================
FMM3D_DEBUG ?= 0
# Submodule location
FMM3D_DIR      ?= external/FMM3D

# Absolute paths to headers and static archive
FMM3D_ROOT     := $(abspath $(FMM3D_DIR))
FMM3D_LIB      := $(FMM3D_ROOT)/local

.PHONY: fmm3d
fmm3d:
ifeq ($(USE_FMM3D),1)
	@echo "==> FMM3D: building via upstream makefile (install)"
	@cd "$(FMM3D_DIR)" && \
	  $(MAKE) install PREFIX=$(abspath $(FMM3D_DIR)/local) DO_DEBUG=$(FMM3D_DEBUG) FAST_KER=OFF
else
	@echo "USE_FMM3D=0 -> skipping FMM3D build"
endif

# Make high-level targets depend on FMM only when enabled
ifeq ($(USE_FMM3D),1)
ALL_DEPS := fmm3d
PY_DEPS  := fmm3d
else
ALL_DEPS :=
PY_DEPS  :=
endif
#=======================================================================

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
			/DUSE_CUDA=${USE_CUDA} /DUSE_MICROMAG=${USE_MICROMAG} /DUSE_FMM3D=${USE_FMM3D}
	else
		FFLAGS = -O3 -fpp -real-size 64 -qopenmp -assume nocc_omp -fpe0 \
			-heap-arrays 1024 -traceback \
			-fp-model=source -fpic -nologo -DUSE_CVODE=${USE_CVODE} \
			-DUSE_MATLAB=${USE_MATLAB} -DUSE_CUDA=${USE_CUDA} \
			-DUSE_MICROMAG=${USE_MICROMAG} -DUSE_FMM3D=${USE_FMM3D}

	endif
else ifeq (${FC}, gfortran)
	FFLAGS = -O3 -fdefault-real-8 -fopenmp -ffree-line-length-512 -cpp -fPIC \
		-DUSE_MICROMAG=0 -DUSE_CVODE=${USE_CVODE} -DUSE_FMM3D=${USE_FMM3D}

endif

PYTHON_MODN = magtensesource
PYTHON_LIBPATH = python/src/magtense/lib

AUX_PATH = source/Aux
NUM_INT_PATH = source/NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = source/TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = source/DemagField/DemagField
MICROMAG_PATH = source/MagTenseMicroMag
FORTRAN_CUDA_PATH = source/MagTenseFortranCuda/cuda
STANDALONE_PATH = source/MagTense_StandAlone/MagTense_StandAlone
FORCEINTEGRATOR_PATH = source/MagneticForceIntegrator/MagneticForceIntegrator

VPATH = ${AUX_PATH}:${NUM_INT_PATH}:${TILE_DEMAG_TENSOR_PATH}:${DEMAG_FIELD_PATH}:\
${MICROMAG_PATH}:${FORTRAN_CUDA_PATH}:${STANDALONE_PATH}:${FORCEINTEGRATOR_PATH}

ifeq ($(OS),Windows_NT)
	CONDA_PATH = $(subst \,/,${CONDA_PREFIX})
	CUDA_ROOT = ${CONDA_PATH}/Library/lib
	MKL = -L${CONDA_PATH}/Library/lib -lmkl_intel_lp64_dll -lmkl_intel_thread_dll \
		-lmkl_core_dll -lmkl_blas95_lp64 -llibiomp5md
	LDFLAGS = '/DEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcmt.lib /LIBPATH:${CONDA_PATH}/Library/lib'
	LIB_SUFFIX = .lib
	PY_MOD_SUFFIX = .pyd
	CVODE_SUFFIX = _static

	ifeq (${FC}, ifx)
		EXTRA_FFLAGS = "${FFLAGS} /assume:underscore /names:lowercase"
		OPT = ${CONDA_PATH}/Library/include \
			-I${CONDA_PATH}/Library/include/intel64/lp64 \
			-I${CONDA_PATH}/opt/compiler/include/intel64
		LIB_OPT = -llibDemagField
	endif
else
 	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -liomp5 -lmkl_blas95_lp64 -lpthread -lm -ldl
	CUDA_ROOT = ${CONDA_PREFIX}/lib
	LDFLAGS =
#	LDFLAGS += '-lstdc++ -liomp5'
	LIB_SUFFIX = .a
	PY_MOD_SUFFIX = .so
	CVODE_SUFFIX =

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
		CP_LIB = cp ${AUX_PATH}/libAux${LIB_SUFFIX} .
		CP_LIB += && cp ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX} .
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
		CP_LIB = cp ${AUX_PATH}/libAux${LIB_SUFFIX} .
		CP_LIB += && cp ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX} .
		CP_LIB += && cp ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX} .
		CP_LIB += && cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
		CP_LIB += && cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
	endif
endif

#===================== FMM3D Integration ===========================
ifeq ($(USE_FMM3D),0)
  FMM3D =
else
  ifeq ($(OS),Windows_NT)
    # Link with full path to the static library (most robust with MS link.exe)
    FMM3D = "$(FMM3D_LIB)/libfmm3d.lib"
  else
    # TODO - should there be some windows-specic handling here?
    #          fx. should it be -llibffm3d instead???

    # tell linker where to find libfmm3d.so
    FMM3D = -L${FMM3D_LIB} -lfmm3d
    # ensure python finds it at runtime
    LDFLAGS += -Wl,-rpath,${FMM3D_LIB}
    # still copy it locally for convenience
    CP_LIB += && cp ${FMM3D_LIB}/libfmm3d${LIB_SUFFIX} .
  endif
endif
#===================================================================

AUX = aux
ifeq ($(OS),Windows_NT)
	LIB_OPT += -llibAux
else
	LIB_OPT += -lAux
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
	CUDA = -L${CUDA_ROOT} -lcublas -lcudart -lcusparse
	CP_LIB += && cp ${FORTRAN_CUDA_PATH}/libCuda${LIB_SUFFIX} .
	ifeq ($(OS),Windows_NT)
		CUDA += -lcuda
		LIB_OPT += -llibCuda
	else
		LIB_OPT += -lCuda
	endif
endif

ifeq ($(USE_CVODE),0)
	CVODE =
else
	CVODE = -L${CVODE_ROOT}/lib -lsundials_core${CVODE_SUFFIX} -lsundials_cvode${CVODE_SUFFIX} \
		-lsundials_fcore_mod -lsundials_fcvode_mod${CVODE_SUFFIX} -lsundials_fnvecserial_mod${CVODE_SUFFIX} \
		-lsundials_fsunmatrixdense_mod${CVODE_SUFFIX} -lsundials_fsunlinsolspgmr_mod${CVODE_SUFFIX}
endif

INCLUDE_OBJ = ${MKFILE_PATH}/${AUX_PATH} \
	-I${MKFILE_PATH}/${NUM_INT_PATH} \
	-I${MKFILE_PATH}/${TILE_DEMAG_TENSOR_PATH} \
	-I${MKFILE_PATH}/${DEMAG_FIELD_PATH} \
	-I${MKFILE_PATH}/${MICROMAG_PATH} \
	-I${MKFILE_PATH}/${FORTRAN_CUDA_PATH} 

PYTHON_MODN_ALL = _${PYTHON_MODN}${PY_MOD_SUFFIX}

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: $(ALL_DEPS) ${AUX} magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} 

standalone: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} standalone

python: $(PY_DEPS) ${AUX} magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${PYTHON_MODN_ALL}

python-win: ${PYTHON_MODN_ALL}

clean:
	cd ${AUX_PATH} && ${MAKE} clean
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	cd $(STANDALONE_PATH) && ${MAKE} clean
	cd $(FORCEINTEGRATOR_PATH) && ${MAKE} clean
	rm -f *${LIB_SUFFIX} *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/*${LIB_SUFFIX} ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build

clean_full:
	cd ${FMM3D_DIR} && ${MAKE} clean
	cd ${AUX_PATH} && ${MAKE} clean
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	cd $(STANDALONE_PATH) && ${MAKE} clean
	cd $(FORCEINTEGRATOR_PATH) && ${MAKE} clean
	rm -f *${LIB_SUFFIX} *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/*${LIB_SUFFIX} ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build

aux:
	cd ${AUX_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"

clean-build:
	rm -f ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build

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

info:
	@echo Using Fortran compiler: $(FC)
	@echo Fortran flags: $(FFLAGS)
	@echo Linker flags: $(LDFLAGS)
	@echo MKL flags: $(MKL)
	@echo CUDA enabled: $(USE_CUDA)
	@echo CVODE enabled: $(USE_CVODE)
	@echo Micromagnetics enabled: $(USE_MICROMAG)
	@echo FMM3D enabled: $(USE_FMM3D)
	@echo MATLAB enabled: $(USE_MATLAB)
	@echo Include paths: -I$(INCLUDE_OBJ)
	@echo Libraries: -L${MKFILE_PATH} ${LIB_OPT}



.PHONY: test
test:
	@command -v fpm >/dev/null 2>&1 || (echo "ERROR: fpm not found. Install with: conda install -c conda-forge fpm" && exit 1)
	@echo "==> Cleaning all fpm test build directories"
	@find tests -mindepth 2 -maxdepth 2 -type d -name build -exec rm -rf {} +
	@echo "==> Running all fpm test harnesses"
	cd tests && FPM_FC=$(FC) FPM_FFLAGS='$(FFLAGS)' bash run_all_fpm_tests.sh




${PYTHON_MODN_ALL}:
	${CP_LIB}
	FC=${FC} FFLAGS=${EXTRA_FFLAGS} LDFLAGS=${LDFLAGS} \
		python -m numpy.f2py -c -m ${PYTHON_MODN} \
		--build-dir ${PYTHON_LIBPATH}/build -I${OPT} -I${INCLUDE_OBJ} \
		-L${MKFILE_PATH} ${LIB_OPT} python/FortranToPythonIO.f90 ${MKL} ${CUDA} ${CVODE} ${FMM3D}
	cp *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/
