
FC = /opt/intel/bin/ifort


FFLAGS = -r8 -O3 -fpe0 -fpp -fp-model source -fpic -I../../NumericalIntegration/NumericalIntegration/ \
	-I../../TileDemagTensor/TileDemagTensor/ -I../../DemagField/DemagField/




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

MagneticForceIntegratorLib.a: $(MagneticForceIntegratorLib)
	ar rc MagneticForceIntegratorLib.a MagForceIO.o MagForceMaxwellStressTensor.o MagneticForce.o


clean:
	rm -f *.o *.x *.mod *.a


MagForceIO.o	:MagForceIO.f90
MagForceMaxwellStressTensor.o	:MagForceMaxwellStressTensor.f90 
MagneticForce.o	:MagneticForce.f90
