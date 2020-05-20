
FC = /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort

NumIntLib = ../../NumericalIntegration/NumericalIntegration/
TileDemagTensorLib = ../../TileDemagTensor/TileDemagTensor/
DemagFieldLib = ../../DemagField/DemagField/


FFLAGS = -r8 -O3 -fpe0 -fp-model source -fpic -I$(NumIntLib) -I$(TileDemagTensorLib) -I$(DemagFieldLib)


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

MagTense.x = StandAloneIO.o MagTense_StandAlone.o

MagTense.x: $(MagTense.x)
	$(FC) $(FFLAGS) $(MagTense.x) -Wl,--start-group $(NumIntLib)libNumericalIntegration.a \
	$(TileDemagTensorLib)libTileDemagTensor.a $(DemagFieldLib)libDemagField.a \
	 -Wl,--end-group -o MagTense.x


clean:
	rm -f *.o *.x *.mod *.a


IO.o	:IO.f90
MagTense_StandAlone.o		:MagTense_StandAlone.f90 StandAloneIO.o





