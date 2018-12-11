
FC = /opt/intel/bin/ifort


FFLAGS = -r8 -O3 -fpe0 -fp-model source -fpic -I../../../NumericalIntegrationLib/NumericalIntegration/NumericalIntegration/ \
	-I../../../SpecialFunctions/SpecialFunctions/ -I../../../TileDemagTensor/TileDemagTensor/TileDemagTensor/


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

MagStatPrism2Lib = $(FC) locate.o MagStatParameters.o MagStatUtil.o MagStat2GetSolution.o IterateMagnetSolution.o

MagStatPrism2lib.a: $(MagStatPrism2Lib)
	ar rc MagStatPrism2lib.a locate.o MagStatParameters.o MagStatUtil.o  MagStat2GetSolution.o IterateMagnetSolution.o


clean:
	rm -f *.o *.x *.mod *.a


locate.o	:locate.for
MagStatParameters.o	:MagStatParameters.f90
MagStatUtil.o	:MagStatUtil.f90 locate.o
MagStat2GetSolution.o	:MagStat2GetSolution.f90
IterateMagnetSolution.o	:IterateMagnetSolution.f90 MagStatParameters.o MagStat2GetSolution.o MagStatUtil.o




