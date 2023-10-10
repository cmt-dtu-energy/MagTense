#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifort
# FC = gfortran

CPP = icx
COMPILE_MICROMAG = 1

USE_MATLAB = 1
MATLAB_INCLUDE = /usr/local/MATLAB/R2021b/extern/include

USE_CUDA = 0

USE_CVODE = 1
CVODE_ROOT= /usr/local/sundials-4.1.0/instdir

MKL_ROOT = /opt/intel/oneapi/mkl/latest 						# Linux - oneapi
# MKL_ROOT = ${CONDA_PREFIX} 									# Linux - conda
# MKL_ROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest" 	# Win - oneapi
# MKL_ROOT = ${CONDA_PREFIX}/Library 							# Win - conda

ifeq (${FC}, ifort)
	ifeq ($(OS),Windows_NT)
		FFLAGS = /nologo /real-size:64 /O3 /assume:nocc_omp /Qopenmp /fpp \
		/fpe:0 /fp:source /fp:precise /libs:static /DUSE_CVODE=${USE_CVODE} \
		/DUSE_CUDA=${USE_CUDA} /DUSE_MATLAB=${USE_MATLAB}
	else
		FFLAGS = -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source \
		-fp-model precise -fpic -libs:static -DUSE_CVODE=${USE_CVODE} \
		-DUSE_CUDA=${USE_CUDA} -DUSE_MATLAB=${USE_MATLAB}
	endif
else ifeq (${FC}, gfortran)
FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp -DUSE_CVODE=${USE_CVODE}
USE_CUDA = 0
COMPILE_MICROMAG= 0
endif

ifeq ($(USE_CUDA),0)
	COMPILE_CUDA =
else
	COMPILE_CUDA = cuda
endif

ifeq ($(COMPILE_MICROMAG),0)
	MICROMAG =
else
	MICROMAG = micromagnetism
endif

ifeq ($(USE_MATLAB),0)
	FORCEINTEGRATOR =
	MATLAB_OPT =
else
	FORCEINTEGRATOR = forceintegrator
	MATLAB_OPT = -I${MATLAB_INCLUDE}
endif
#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} # standalone 

magnetostatic:
	cd source/NumericalIntegration/NumericalIntegration && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)' USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) MATLAB_OPT=$(MATLAB_OPT)
	cd source/TileDemagTensor/TileDemagTensor && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd source/DemagField/DemagField && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)' USE_MATLAB=$(USE_MATLAB) MATLAB_OPT=$(MATLAB_OPT)

micromagnetism:
	cd source/MagTenseMicroMag && $(MAKE) FFLAGS='$(FFLAGS)' USE_MATLAB=$(USE_MATLAB) MATLAB_OPT=$(MATLAB_OPT) MKL_ROOT=$(MKL_ROOT)

cuda:
	cd source/MagTenseFortranCuda/cuda && $(MAKE) CPP=$(CPP)

forceintegrator:
	cd source/MagneticForceIntegrator/MagneticForceIntegrator && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)' MATLAB_OPT=$(MATLAB_OPT)

standalone:
	cd source/MagTense_StandAlone/MagTense_StandAlone && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	mkdir build
	cp source/MagTense_StandAlone/MagTense_StandAlone/MagTense.x build/MagTense.x

clean:
	cd source/NumericalIntegration/NumericalIntegration && make clean
	cd source/TileDemagTensor/TileDemagTensor && make clean
	cd source/DemagField/DemagField && make clean
	cd source/MagTenseMicroMag && make clean
	cd source/MagTenseFortranCuda/cuda && make clean
	cd source/MagTense_StandAlone/MagTense_StandAlone && make clean
	cd source/MagneticForceIntegrator/MagneticForceIntegrator && make clean
	rm -r build/
