#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifx
# FC = gfortran
CPP = icx

USE_CUDA = 1
USE_CVODE = 0
USE_MATLAB = 1
COMPILE_MICROMAG = 1

CVODE_ROOT= /usr/local/sundials-4.1.0/instdir
# MATLAB_INCLUDE = /usr/local/MATLAB/R2023b/extern/include			# Linux
MATLAB_INCLUDE = "C:\Program Files\MATLAB\R2023b\extern\include"	# Win
MKL_ROOT = ${CONDA_PREFIX} 											# Linux - mkl
# MKL_ROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest" 		# Win - mkl (oneapi)
# MKL_ROOT = ${CONDA_PREFIX}/Library 								# Win - mkl (conda)

ifeq (${FC}, ifort)
	ifeq ($(OS),Windows_NT)
		FFLAGS = /O3 /fpp /real-size:64 /assume:nocc_omp /Qopenmp \
		/fpe:0 /fp:source /libs:static /DUSE_CVODE=${USE_CVODE} \
		/DUSE_CUDA=${USE_CUDA} /DUSE_MATLAB=${USE_MATLAB}
	else
		FFLAGS = -O3 -fpp -real-size 64 -assume nocc_omp -qopenmp \
		-fpe0 -fp-model=source -fpic -DUSE_CVODE=${USE_CVODE} \
		-DUSE_CUDA=${USE_CUDA} -DUSE_MATLAB=${USE_MATLAB}
	endif
else ifeq (${FC}, gfortran)
FFLAGS = -O3 -fdefault-real-8 -fopenmp -ffree-line-length-512 -cpp -fPIC -DUSE_CVODE=${USE_CVODE}
USE_CUDA = 0
COMPILE_MICROMAG= 0
endif

ifeq ($(USE_CUDA),0)
	COMPILE_CUDA =
else
	COMPILE_CUDA = cuda
endif

ifeq ($(USE_MATLAB),0)
	FORCEINTEGRATOR =
else
	FORCEINTEGRATOR = forceintegrator
endif

ifeq ($(COMPILE_MICROMAG),0)
	MICROMAG =
else
	MICROMAG = micromagnetism
endif
#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} # standalone

magnetostatic:
	cd source/NumericalIntegration/NumericalIntegration && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd source/TileDemagTensor/TileDemagTensor && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)
	cd source/DemagField/DemagField && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE)

micromagnetism:
	cd source/MagTenseMicroMag && $(MAKE) FFLAGS="$(FFLAGS)" USE_CVODE=$(USE_CVODE) CVODE_ROOT=$(CVODE_ROOT) USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE=$(MATLAB_INCLUDE) MKL_ROOT=$(MKL_ROOT)

cuda:
	cd source/MagTenseFortranCuda/cuda && $(MAKE) CPP=$(CPP)

forceintegrator:
	cd source/MagneticForceIntegrator/MagneticForceIntegrator && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" MATLAB_INCLUDE=$(MATLAB_INCLUDE)

standalone:
	cd source/MagTense_StandAlone/MagTense_StandAlone && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)"
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
