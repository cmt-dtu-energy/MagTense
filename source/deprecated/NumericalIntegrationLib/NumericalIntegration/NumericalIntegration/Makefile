
FC = /opt/intel/bin/ifort


FFLAGS = -r8 -O3 -fpe0 -fp-model source -fpic


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

NumericalIntegrationLib = $(FC) IntegrationDataTypes.o quadpack.o NumericalIntegrals.o

NumericalIntegrationlib.a: $(NumericalIntegrationLib)
	ar rc NumericalIntegrationlib.a IntegrationDataTypes.o quadpack.o NumericalIntegrals.o


clean:
	rm -f *.o *.x *.mod *.a


quadpack.o	:quadpack.f90 IntegrationDataTypes.o
IntegrationDataTypes.o	:IntegrationDataTypes.f90 
NumericalIntegrals.o	:NumericalIntegrals.f90 quadpack.o IntegrationDataTypes.o

