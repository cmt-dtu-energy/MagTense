#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: magnetostatic # micromag standalone forceintegrator

magnetostatic:
	cd source/NumericalIntegration/NumericalIntegration && $(MAKE)
	cd source/TileDemagTensor/TileDemagTensor && $(MAKE)
	cd source/DemagField/DemagField && $(MAKE)

micromag:
	cd source/MagTenseMicroMag && $(MAKE)

forceintegrator:
	cd source/MagneticForceIntegrator/MagneticForceIntegrator && $(MAKE)

standalone:
	cd source/MagTense_StandAlone/MagTense_StandAlone && $(MAKE)
	mkdir build
	cp source/MagTense_StandAlone/MagTense_StandAlone/MagTense.x build/MagTense.x

clean:
	cd source/NumericalIntegration/NumericalIntegration && make clean
	cd source/TileDemagTensor/TileDemagTensor && make clean
	cd source/DemagField/DemagField && make clean
	cd source/MagTenseMicroMag && make clean
	cd source/MagTense_StandAlone/MagTense_StandAlone && make clean
	cd source/MagneticForceIntegrator/MagneticForceIntegrator && make clean
