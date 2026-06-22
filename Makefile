#=======================================================================
#                       compiler names and flags
#=======================================================================
USE_CUDA = 1
PY_VERSION ?= 314
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

# Location where Miniconda will be installed
CONDA_DIR := $(HOME)/miniconda3
CONDA_BIN := $(CONDA_DIR)/bin/conda
#=======================================================================
#                    FMM3D integration (upstream Makefile)
#=======================================================================
FMM3D_DEBUG ?= 0
# Submodule location
FMM3D_DIR      ?= external/FMM3D

# Absolute paths to headers and static archive
FMM3D_ROOT     := $(abspath $(FMM3D_DIR))
FMM3D_LIB      := $(FMM3D_ROOT)/local

ifeq ($(OS),Windows_NT)
  SEP = ;
else
  SEP = &&
endif

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

# OS-specific shell
ifeq ($(OS),Windows_NT)
  SHELL := powershell.exe
  .SHELLFLAGS := -NoProfile -ExecutionPolicy Bypass -Command
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

AUXMT_PATH = source/AuxMT
NUM_INT_PATH = source/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = source/TileDemagTensor
DEMAG_FIELD_PATH = source/DemagField
MICROMAG_PATH = source/MagTenseMicroMag
FORTRAN_CUDA_PATH = source/MagTenseFortranCuda/cuda
STANDALONE_PATH = source/MagTense_StandAlone
FORCEINTEGRATOR_PATH = source/MagneticForceIntegrator

VPATH = ${AUXMT_PATH}:${NUM_INT_PATH}:${TILE_DEMAG_TENSOR_PATH}:${DEMAG_FIELD_PATH}:\
${MICROMAG_PATH}:${FORTRAN_CUDA_PATH}:${STANDALONE_PATH}:${FORCEINTEGRATOR_PATH}

ifeq ($(OS),Windows_NT)
	CONDA_PATH = $(subst \,/,${CONDA_PREFIX})
	CUDA_ROOT = ${CONDA_PATH}/Library/lib
	MKL = -L${CONDA_PATH}/Library/lib -lmkl_intel_lp64_dll -lmkl_intel_thread_dll \
		-lmkl_core_dll -lmkl_blas95_lp64 -llibiomp5md
		
	AUXMT_PATH_ROOT     := $(abspath $(AUXMT_PATH))
	
	ifeq ($(USE_FMM3D),0)	
		LDFLAGS = '/DEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcmt.lib /LIBPATH:${CONDA_PATH}/Library/lib'
	else
		LDFLAGS = '/DEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcmt.lib /LIBPATH:${CONDA_PATH}/Library/lib /LIBPATH:${FMM3D_LIB}'
	endif
	
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


AUXMT = auxmt
ifeq ($(OS),Windows_NT)
	LIB_OPT += -llibAuxMT
else
	LIB_OPT += -lAuxMT
endif

ifeq ($(USE_MICROMAG),0)
	MICROMAG =
	ifeq ($(OS),Windows_NT)
		CP_LIB = cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
	else
		CP_LIB = cp ${AUXMT_PATH}/libAuxMT${LIB_SUFFIX} .
		CP_LIB += && cp ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX} .
		CP_LIB += && cp ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX} .
		CP_LIB += && cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
	endif
else
	MICROMAG = micromagnetism
	ifeq ($(OS),Windows_NT)
		LIB_OPT = -llibAuxMT
		CP_LIB = cp ${AUXMT_PATH}/libAuxMT${LIB_SUFFIX} .
		
		LIB_OPT += -llibMagTenseMicroMag
		#CP_LIB = cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
		CP_LIB += && cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
		
	else
		LIB_OPT += -lMagTenseMicroMag
		CP_LIB = cp ${AUXMT_PATH}/libAuxMT${LIB_SUFFIX} .
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
    # 1. Use the standard linker flag (-lfmm3d looks for fmm3d.lib)
    FMM3D = -L"${FMM3D_LIB}" -lfmm3d
    
    # 2. Use 'cp' (since your shell is sh) and RENAME the file during copy
    # This creates fmm3d.lib in your root so the linker finds it easily
    CP_LIB += && cp "${FMM3D_LIB}/libfmm3d.lib" ./fmm3d.lib
  else
    FMM3D = -L${FMM3D_LIB} -lfmm3d
    LDFLAGS += -Wl,-rpath,${FMM3D_LIB}
    CP_LIB += && cp ${FMM3D_LIB}/libfmm3d${LIB_SUFFIX} .
  endif
endif
#===================================================================



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

INCLUDE_OBJ = ${MKFILE_PATH}/${AUXMT_PATH} \
	-I${MKFILE_PATH}/${NUM_INT_PATH} \
	-I${MKFILE_PATH}/${TILE_DEMAG_TENSOR_PATH} \
	-I${MKFILE_PATH}/${DEMAG_FIELD_PATH} \
	-I${MKFILE_PATH}/${MICROMAG_PATH} \
	-I${MKFILE_PATH}/${FORTRAN_CUDA_PATH} 

PYTHON_MODN_ALL = _${PYTHON_MODN}${PY_MOD_SUFFIX}

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean install-miniconda build-env build-cvode ${MICROMAG_PATH} test standalone python_ python python-win ${PYTHON_MODN_ALL} rm-conda rm-env python-interface-win 

all: $(ALL_DEPS) ${AUXMT} magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} 

standalone: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR}

python: $(PY_DEPS) ${AUXMT} magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${PYTHON_MODN_ALL}

python-win: ${PYTHON_MODN_ALL}

define clean-subdirs
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	cd $(STANDALONE_PATH) && ${MAKE} clean
	cd $(FORCEINTEGRATOR_PATH) && ${MAKE} clean
endef

clean:
	$(clean-subdirs)
	rm -f *${LIB_SUFFIX} *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/*${LIB_SUFFIX} ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build
	rm -rf cvode*

clean_full:
	cd ${FMM3D_DIR} && ${MAKE} clean
	cd ${AUXMT_PATH} && ${MAKE} clean
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	cd $(STANDALONE_PATH) && ${MAKE} clean
	cd $(FORCEINTEGRATOR_PATH) && ${MAKE} clean
	rm -f *${LIB_SUFFIX} *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/*${LIB_SUFFIX} ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build

auxmt:
	cd ${AUXMT_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"

clean-build:
	rm -f ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build

magnetostatic:
	cd ${NUM_INT_PATH} $(SEP) $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd ${TILE_DEMAG_TENSOR_PATH} $(SEP) $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd ${DEMAG_FIELD_PATH}  $(SEP) $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)

micromagnetism:
	cd ${MICROMAG_PATH} $(SEP) $(MAKE) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)

cuda:
	cd ${FORTRAN_CUDA_PATH} $(SEP) $(MAKE) CPP=$(CPP)

forceintegrator:
	cd $(FORCEINTEGRATOR_PATH) $(SEP) $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' MATLAB_INCLUDE=$(MATLAB_INCLUDE)

standalone:
	cd $(STANDALONE_PATH) $(SEP) $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}'
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

# Rule that installs Miniconda only if "conda" is not found
install-miniconda:
	@if [ -x "$(CONDA_BIN)" ]; then \
			echo "Conda already installed at $(CONDA_DIR)."; \
	else \
			echo "Conda not found at $(CONDA_DIR). Installing Miniconda..."; \
			curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh; \
			bash miniconda.sh -b -p $(CONDA_DIR) && rm miniconda.sh; \
			echo "Miniconda installed at $(CONDA_DIR)."; \
	fi

rm-conda:
	rm -rf $(CONDA_DIR)

build-env: install-miniconda
# Check if the "magtense-env" environment already exists before creating it
	@if ! conda env list  | grep -q "magtense-env"; then \
		$(CONDA_BIN) env create -n magtense-env -f ${MKFILE_PATH}/python/.build/env-$(PY_VERSION)-linux.yml; \
	else \
		echo "magtense-env environment already exists."; \
	fi


rm-env:
	$(CONDA_BIN) env remove -n magtense-env -y

CMAKE = $(CONDA_BIN) run -n magtense-env -- cmake
IFX = $(CONDA_BIN) run -n magtense-env which ifx
ICX = $(CONDA_BIN) run -n magtense-env which icx


define run-cmake-cvode
	$(CMAKE) \
	-B ${CVODE_ROOT}/build \
	-S ${CVODE_ROOT}/src \
	-D CMAKE_BUILD_TYPE=Release \
	-D BUILD_ARKODE=OFF \
	-D BUILD_CVODE=ON \
	-D BUILD_CVODES=OFF \
	-D BUILD_IDA=OFF \
	-D BUILD_IDAS=OFF \
	-D BUILD_KINSOL=OFF \
	-D BUILD_SHARED_LIBS=OFF \
	-D BUILD_STATIC_LIBS=ON \
	-D CMAKE_INSTALL_PREFIX=${CVODE_ROOT} \
	-D EXAMPLES_INSTALL_PATH=${CVODE_ROOT}/examples \
	-D CMAKE_C_COMPILER=$$($(ICX)) \
	-D CMAKE_Fortran_COMPILER=$$($(IFX)) \
	-D BUILD_FORTRAN_MODULE_INTERFACE=ON \
	-D ENABLE_OPENMP=ON
	$(CMAKE) --build ${CVODE_ROOT}/build --config Release --verbose
	$(CMAKE) --install ${CVODE_ROOT}/build --verbose
endef

build-cvode: build-env
	wget https://github.com/LLNL/sundials/releases/download/v7.4.0/cvode-7.4.0.tar.gz
	tar -xf cvode-7.4.0.tar.gz
	rm -rf ${CVODE_ROOT}
	mkdir -p ${CVODE_ROOT}
	mv ${MKFILE_PATH}/cvode-7.4.0 ${CVODE_ROOT}/src
	$(run-cmake-cvode)

ENV_NAME := magtense-env

python-interface: build-env
	@if [ ! -d "${CVODE_ROOT}/build" ] || [ -z "$$(ls -A ${CVODE_ROOT}/build)" ]; then \
		$(MAKE) build-cvode; \
	fi
	cp python/.build/requirements-py3-dev.txt python/requirements.txt
	
	@if [ "$$CONDA_DEFAULT_ENV" = "$(ENV_NAME)" ]; then \
		$(MAKE) python USE_CUDA=$(USE_CUDA) USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=0 && \
		python -m pip install -e ./python; \
	else \
		$(CONDA_BIN) run -n magtense-env $(MAKE) python USE_CUDA=$(USE_CUDA) USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=0 && \
		$(CONDA_BIN) run -n magtense-env -- python -m pip install -e ./python; \
	fi

pytest:
	$(CONDA_BIN) run -n magtense-env -- pytest

python-interface-win: 
	make magnetostatic micromagnetism USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=0
