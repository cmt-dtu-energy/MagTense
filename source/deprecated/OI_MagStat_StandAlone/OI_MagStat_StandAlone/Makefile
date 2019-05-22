
FC = /opt/intel/bin/ifort

NumIntLib = ../../NumericalIntegrationLib/NumericalIntegration/NumericalIntegration/
SpecFuncLib = ../../SpecialFunctions/SpecialFunctions/
DemagTensorLib = ../../TileDemagTensor/TileDemagTensor/TileDemagTensor/
MagStatLib = ../../MagStatPrisms2/MagStatVersion2/MagStatVersion2/


FFLAGS = -r8 -O3 -fpe0 -fp-model source -fpic -I$(NumIntLib) -I$(SpecFuncLib) -I$(DemagTensorLib) -I$(MagStatLib)


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

MagStatPrism.x = IO.o OI_MagStat_StandAlone.o

MagStatPrism.x: $(MagStatPrism.x)
	$(FC) $(FFLAGS) $(MagStatPrism.x) -Wl,--start-group $(NumIntLib)NumericalIntegrationlib.a \
	$(SpecFuncLib)special_functions_lib.a $(DemagTensorLib)TileDemagTensorlib.a $(MagStatLib)MagStatPrism2lib.a \
	 -Wl,--end-group -o MagStatPrism.x


clean:
	rm -f *.o *.x *.mod *.a


IO.o	:IO.f90
OI_MagStat_StandAlone.o		:OI_MagStat_StandAlone.f90 IO.o





