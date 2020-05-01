
FC = /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort


FFLAGS = -r8 -O3 -fpp -fpe0 -fp-model source -fpic -I../NumericalIntegration/NumericalIntegration/ \
	-I../DemagField/DemagField/ -I/apps/external/matlab/2019/b/extern/include/


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

MagTenseMicroMagLib = $(FC) LandauLifshitzEquationSolver.o MagTenseMicroMagIO.o MicroMagParameters.o FortranCuda.o LLODE_Debug.o

MagTenseMicroMagLib.a: $(MagTenseMicroMagLib)
	ar rc DemagFieldLib.a DemagFieldGetSolution.o MagTenseMicroMagIO.o MicroMagParameters.o FortranCuda.o LLODE_Debug.o


clean:
	rm -f *.o *.x *.mod *.a


MagTenseMicroMagIO.o	:MagTenseMicroMagIO.f90
MicroMagParameters.o	:MicroMagParameters.f90
FortranCuda.o	:FortranCuda.f90
LLODE_Debug.o	:LLODE_Debug.f90 MicroMagParameters.o
LandauLifshitzEquationSolver.o	:LandauLifshitzEquationSolver.f90 LLODE_Debug.o FortranCuda.o MicroMagParameters.o MagTenseMicroMagIO.o



