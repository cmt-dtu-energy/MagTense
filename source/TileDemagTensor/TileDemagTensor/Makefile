
FC = /apps/external/intel/2019/compilers_and_libraries_2019.4.243/linux/bin/intel64/ifort


FFLAGS = -r8 -O3 -fpe0 -fp-model source -fpic -I../../NumericalIntegration/NumericalIntegration/ \


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

TileDemagTensorLib = $(FC) SpecialFunctions.o TileTensorHelperFunctions.o TileTriangle.o TileCircPieceTensor.o TileCylPieceTensor.o TileRectangularPrismTensor.o TilePlanarCoilTensor.o  TileNComponents.o

TileDemagTensorlib.a: $(TileDemagTensorLib)
	ar rc TileDemagTensorlib.a SpecialFunctions.o TileTensorHelperFunctions.o TileTriangle.o TileCircPieceTensor.o TileCylPieceTensor.o TileRectangularPrismTensor.o TilePlanarCoilTensor.o  TileNComponents.o


clean:
	rm -f *.o *.x *.mod *.a

SpecialFunctions.o	:SpecialFunctions.f90
TileTensorHelperFunctions.o	:TileTensorHelperFunctions.f90
TileCircPieceTensor.o	:TileCircPieceTensor.f90 TileTensorHelperFunctions.o
TileCylPieceTensor.o   :TileCylPieceTensor.f90 TileTensorHelperFunctions.o
TileRectangularPrismTensor.o   :TileRectangularPrismTensor.f90 TileTensorHelperFunctions.o
TilePlanarCoilTensor.o   :TilePlanarCoilTensor.f90 TileTensorHelperFunctions.o
TileTriangle.o   :TileTriangle.f90 TileTensorHelperFunctions.o
TileNComponents.o	:TileNComponents.f90 TileTensorHelperFunctions.o TileTriangle.o TileCircPieceTensor.o TileCylPieceTensor.o TileRectangularPrismTensor.o TilePlanarCoilTensor.o




