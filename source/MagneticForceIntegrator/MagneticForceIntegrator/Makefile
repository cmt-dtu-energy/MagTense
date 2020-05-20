
FC = /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort


FFLAGS = -r8 -O3 -fpe0 -fpp -fp-model source -fpic -static-intel -libs:static -I../../NumericalIntegration/NumericalIntegration/ \
	-I../../TileDemagTensor/TileDemagTensor/ -I../../DemagField/DemagField/  -I/apps/external/matlab/2019/b/extern/include/ 




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

MagneticForceIntegratorLib = $(FC) MagForceIO.o MagForceMaxwellStressTensor.o MagneticForce.o

libMagneticForceIntegrator.a: $(MagneticForceIntegratorLib)
	/apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/xiar rc libMagneticForceIntegrator.a MagForceIO.o MagForceMaxwellStressTensor.o MagneticForce.o


clean:
	rm -f *.o *.x *.mod *.a


MagForceIO.o	:MagForceIO.f90
MagForceMaxwellStressTensor.o	:MagForceMaxwellStressTensor.f90 
MagneticForce.o	:MagneticForce.f90
