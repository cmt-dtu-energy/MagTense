
FC = /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort


FFLAGS = -r8 -O3 -fpp -fpe0 -fp-model source -fpic


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

NumericalIntegrationLib = $(FC) IntegrationDataTypes.o quadpack.o NumericalIntegrals.o ODE_Solvers.o ODE_MagTense.o RKSuitef90.o util.o

NumericalIntegrationlib.a: $(NumericalIntegrationLib)
	ar rc NumericalIntegrationlib.a IntegrationDataTypes.o quadpack.o NumericalIntegrals.o ODE_Solvers.o ODE_MagTense.o RKSuitef90.o util.o


clean:
	rm -f *.o *.x *.mod *.a


quadpack.o	:quadpack.f90 IntegrationDataTypes.o
IntegrationDataTypes.o	:IntegrationDataTypes.f90 
NumericalIntegrals.o	:NumericalIntegrals.f90 quadpack.o IntegrationDataTypes.o
RKSuitef90.o	:RKSuitef90.f90
ODE_Solvers.o	:ODE_Solvers.f90 RKSuitef90.o
util.o	:util.f90
