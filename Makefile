#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran

USE_MATLAB = 0
USE_CUDA = 1
USE_CVODE = 0

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
endif
#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: magnetostatic # micromag standalone forceintegrator

magnetostatic:
	cd source/NumericalIntegration/NumericalIntegration && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd source/TileDemagTensor/TileDemagTensor && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd source/DemagField/DemagField && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'

micromag:
	cd source/MagTenseMicroMag && $(MAKE) FFLAGS='$(FFLAGS)' USE_MATLAB=$(USE_MATLAB)

forceintegrator:
	cd source/MagneticForceIntegrator/MagneticForceIntegrator && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'

standalone:
	cd source/MagTense_StandAlone/MagTense_StandAlone && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	mkdir build
	cp source/MagTense_StandAlone/MagTense_StandAlone/MagTense.x build/MagTense.x

clean:
	cd source/NumericalIntegration/NumericalIntegration && make clean
	cd source/TileDemagTensor/TileDemagTensor && make clean
	cd source/DemagField/DemagField && make clean
	cd source/MagTenseMicroMag && make clean
	cd source/MagTense_StandAlone/MagTense_StandAlone && make clean
	cd source/MagneticForceIntegrator/MagneticForceIntegrator && make clean
	rm -r build/
