
FC = /opt/intel/bin/ifort


FFLAGS = -r8 -O3 -fpe0 -fp-model source -fpic -I../../NumericalIntegration/NumericalIntegration/ \


#FC = /usr/bin/gfortran

#FFLAGS = -O3 -fopenmp -fdefault-real-8 -ffixed-line-length-none -fPIC




.KEEP_STATE:
.SUFFIXES:
.SUFFIXES: .for .f90 .F90 .f .o
.f.o:
	$(FC) -c $(FFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $<
.for.o:
	$(FC) -c $(FFLAGS) $<

all: 
	cd source/NumericalIntegration/NumericalIntegration&& $(MAKE) 
	cd source/TileDemagTensor/TileDemagTensor && $(MAKE)
	cd source/DemagField/DemagField && $(MAKE) 
	cd source/MagTense_StandAlone/MagTense_StandAlone && $(MAKE) 
	cp source/MagTense_StandAlone/MagTense_StandAlone/MagTense.x executable/MagTense.x 

clean:
	rm -f *.o *.x *.mod *.a

