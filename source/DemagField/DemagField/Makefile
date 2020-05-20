
FC = /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort


FFLAGS = -r8 -O3 -fpp -fpe0 -fp-model source -fpic -static-intel -libs:static -I../../NumericalIntegration/NumericalIntegration/ \
	-I../../TileDemagTensor/TileDemagTensor/ -I/apps/external/matlab/2019/b/extern/include/ 


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

libDemagField.a: $(DemagFieldLib)
	/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/xiar rc libDemagField.a DemagFieldGetSolution.o MagParameters.o MagTileIO.o  spline.o IterateMagnetSolution.o


clean:
	rm -f *.o *.x *.mod *.a


DemagFieldGetSolution.o	:DemagFieldGetSolution.f90
MagParameters.o	:MagParameters.f90
MagTileIO.o	:MagTileIO.f90
spline.o	:spline.f90
IterateMagnetSolution.o	:IterateMagnetSolution.f90 MagParameters.o DemagFieldGetSolution.o spline.o MagTileIO.o




