
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

SpecialFunctionsLib = $(FC) SpecialFunctions.o

special_functions_lib.a: $(SpecialFunctionsLib)
	ar rc special_functions_lib.a SpecialFunctions.o


clean:
	rm -f *.o *.x *.mod *.a


SpecialFunctions.o	:SpecialFunctions.f90

