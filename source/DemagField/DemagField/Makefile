
FC = /opt/intel/bin/ifort


FFLAGS = -r8 -O3 -fpp -fpe0 -fp-model source -fpic -I../../NumericalIntegration/NumericalIntegration/ \
	-I../../TileDemagTensor/TileDemagTensor/


#FC = /usr/bin/gfortran

#FFLAGS = -O3 -fpp -fopenmp -fdefault-real-8 -ffixed-line-length-none -fPIC




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

DemagFieldLib = $(FC) DemagFieldGetSolution.o MagParameters.o MagTileIO.o spline.o IterateMagnetSolution.o

DemagFieldLib.a: $(DemagFieldLib)
	ar rc DemagFieldLib.a DemagFieldGetSolution.o MagParameters.o MagTileIO.o  spline.o IterateMagnetSolution.o


clean:
	rm -f *.o *.x *.mod *.a


DemagFieldGetSolution.o	:DemagFieldGetSolution.f90
MagParameters.o	:MagParameters.f90
MagTileIO.o	:MagTileIO.f90
spline.o	:spline.f90
IterateMagnetSolution.o	:IterateMagnetSolution.f90 MagParameters.o DemagFieldGetSolution.o spline.o MagTileIO.o




