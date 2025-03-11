#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifx
# FC = gfortran
CPP = icx

USE_CUDA = 1
USE_CVODE = 0
USE_MATLAB = 1
USE_MICROMAG = 1

CVODE_ROOT = /usr/local/sundials-4.1.0/instdir
# MATLAB_INCLUDE = /usr/local/MATLAB/R2023b/extern/include			# Linux
MATLAB_INCLUDE = "C:\Program Files\MATLAB\R2023b\extern\include"	# Win
MKL_ROOT = ${CONDA_PREFIX} 											# Linux - mkl
# MKL_ROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest" 		# Win - mkl (oneapi)
# MKL_ROOT = ${CONDA_PREFIX}/Library 								# Win - mkl (conda)

ifeq (${FC}, ifx)
	ifeq ($(OS),Windows_NT)
		FFLAGS = /O3 /fpp /real-size:64 /Qopenmp /assume:nocc_omp \
		/fpe:0 /fp:source /nologo /libs:static /DUSE_CVODE=${USE_CVODE} \
		/DUSE_MATLAB=${USE_MATLAB} /DUSE_MICROMAG=${USE_MICROMAG}

		LIB_SUFFIX = .lib
	else
		FFLAGS = -O3 -fpp -static-intel -real-size 64 -assume nocc_omp -qopenmp \
		-fpe0 -fp-model=source -fPIC -nologo -DUSE_CVODE=${USE_CVODE} \
		-DUSE_MATLAB=${USE_MATLAB} -DUSE_MICROMAG=${USE_MICROMAG}

		LIB_SUFFIX = .a
	endif
else ifeq (${FC}, gfortran)
	FFLAGS = -O3 -fdefault-real-8 -fopenmp -ffree-line-length-512 -cpp -fPIC \
		-DUSE_CVODE=${USE_CVODE} -DUSE_MICROMAG=0 -DUSE_CVODE=${USE_CVODE}
endif

NUM_INT_PATH = source/NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = source/TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = source/DemagField/DemagField
MICROMAG_PATH = source/MagTenseMicroMag
FORTRAN_CUDA_PATH = source/MagTenseFortranCuda/cuda
STANDALONE_PATH = source/MagTense_StandAlone/MagTense_StandAlone
FORCEINTEGRATOR_PATH = source/MagneticForceIntegrator/MagneticForceIntegrator

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
else
	COMPILE_CUDA = cuda
	CP_LIB += && cp ${FORTRAN_CUDA_PATH}/libCuda${LIB_SUFFIX} .
endif
#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} # standalone

magnetostatic:
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd ${DEMAG_FIELD_PATH}  && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)

micromagnetism:
	cd ${MICROMAG_PATH} && $(MAKE) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)

cuda:
	cd ${FORTRAN_CUDA_PATH} && $(MAKE) CPP=$(CPP)

forceintegrator:
	cd $(FORCEINTEGRATOR_PATH) && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" MATLAB_INCLUDE=$(MATLAB_INCLUDE)

standalone:
	cd $(STANDALONE_PATH) && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)"
	mkdir build
	cp $(STANDALONE_PATH)/MagTense.x build/MagTense.x

clean:
	cd ${NUM_INT_PATH} && make clean
	cd ${TILE_DEMAG_TENSOR_PATH} && make clean
	cd ${DEMAG_FIELD_PATH}  && make clean
	cd ${MICROMAG_PATH} && make clean
	cd ${FORTRAN_CUDA_PATH} && make clean
	cd $(STANDALONE_PATH) && make clean
	cd $(FORCEINTEGRATOR_PATH) && make clean
	rm *${LIB_SUFFIX}
