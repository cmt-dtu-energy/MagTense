#==========================
# Compiler names and flags
#==========================
USE_CUDA = 1
USE_MICROMAG = 1

ifeq ($(OS),Windows_NT)
	F90 = ifort
else
	ifeq (${UNAME}, Darwin)
  		F90 = gfortran
	else
  		F90 = ifort
	endif
endif

ifeq (${F90}, gfortran)
	USE_CUDA = 0
	USE_MKL = 0
	USE_MICROMAG = 0
endif

PYTHON_MODN = magtensesource

NUM_INT_PATH = ../../../../source/NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../../../../source/TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = ../../../../source/DemagField/DemagField
MICROMAG_PATH = ../../../../source/MagTenseMicroMag
FORTRAN_CUDA_PATH = ../../../../source/MagTenseFortranCuda/cuda

VPATH = ${NUM_INT_PATH}:${TILE_DEMAG_TENSOR_PATH}:${DEMAG_FIELD_PATH}:\
${MICROMAG_PATH}:${FORTRAN_CUDA_PATH}

NUM_INT_LIB_NAME = libNumericalIntegration
TILE_DEMAG_TENSOR_LIB_NAME = libTileDemagTensor
DEMAG_FIELD_LIB_NAME = libDemagField
MICROMAG_LIB_NAME = libMagTenseMicroMag

ifeq ($(OS),Windows_NT)
	MKL_ROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
	CUDA = -L${CONDA_PREFIX}/Lib/x64 -lcublas -lcudart -lcuda -lcusparse
	MKL = -L${MKL_ROOT}/lib/intel64 -lmkl_rt -lmkl_blas95_ilp64
	MKL_OMP = -L${MKL_ROOT}/windows/compiler/lib/intel64_win -llibiomp5md

	ifeq (${F90}, ifort)
		F90FLAGS = /nologo /O3 /assume:nocc_omp /Qopenmp /real-size:64 \
		/fp:precise /libs:static /threads /fpp /DUSE_MICROMAG=${USE_MICROMAG}
		FCOMP = intelvem
		EXTRA_OPT = -DNO_APPEND_FORTRAN -DUPPERCASE_FORTRAN --no-lower \
		--opt='${F90FLAGS} /assume:nounderscore /names:uppercase'
		OPT = ${CONDA_PREFIX}/Library/include -I${MKL_ROOT}/include \
		-I${CONDA_PREFIX}/Library/include/intel64/ilp64 -I${NUM_INT_PATH} \
		-I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} \
		-I${MICROMAG_PATH} -I${FORTRAN_CUDA_PATH}
		LIB_OPT = -I:${OPT} -L. -llibsrc
	else ifeq (${F90}, gfortran)
		F90FLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp -DUSE_MICROMAG=0
		FCOMP = gnu95 --compiler=mingw32
		EXTRA_OPT = --opt='${F90FLAGS}'
		LIB_OPT = -L. -llibsrc
	endif
		
	LIBTOOL = LIB $? /OUT:$@
	STATIC_LIB = libsrc.lib
	NUM_INT_LIB = ${NUM_INT_PATH}/${NUM_INT_LIB_NAME}.lib
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/${TILE_DEMAG_TENSOR_LIB_NAME}.lib
	DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/${DEMAG_FIELD_LIB_NAME}.lib
	MICROMAG_LIB = ${MICROMAG_PATH}/${MICROMAG_LIB_NAME}.lib
	CLEAN_FILES = -del libsrc.lib *.o *.mod *.pyd
	CLEAN_FOLDER = -rmdir /s /q Release 2>NUL || VER>NUL
	CLEAN_FOLDER_WIN = -for /f %%i in ('dir /a:d /b src.win-*') do rd /s /q %%i
	PYTHON_MODN_ALL = _${PYTHON_MODN}.dll

else
	ifeq (${F90}, ifort)
		F90FLAGS = -fPIC -nologo -O3 -assume nocc_omp -qopenmp -real-size 64 \
		-fp-model precise -threads -fpp -DUSE_MICROMAG=${USE_MICROMAG}
		FCOMP = intelem
		OPT = ${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/intel64/ilp64 \
		-I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} \
		-I${MICROMAG_PATH} -I${FORTRAN_CUDA_PATH}
		LIB_OPT = -I${OPT} -L. -lsrc
	else ifeq (${F90}, gfortran)
		F90FLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp -DUSE_MICROMAG=0
		FCOMP = gnu95
		LIB_OPT = -L. -lsrc
	endif

	EXTRA_OPT = --opt='${F90FLAGS}'
	CUDA = -L${CONDA_PREFIX}/lib -lcublas -lcudart -lcuda -lcusparse
	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -lpthread -lm -ldl -liomp5 -lmkl_blas95_ilp64

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	
	STATIC_LIB = libsrc.a
	NUM_INT_LIB = ${NUM_INT_PATH}/${NUM_INT_LIB_NAME}.a
	TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/${TILE_DEMAG_TENSOR_LIB_NAME}.a
	DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/${DEMAG_FIELD_LIB_NAME}.a
	MICROMAG_LIB = ${MICROMAG_PATH}/${MICROMAG_LIB_NAME}.a
	CLEAN_FILES = -rm libsrc.a *.o *.mod *.so 
	CLEAN_FOLDER = -rm -r src.*
	CLEAN_FOLDER_WIN = -true
	PYTHON_MODN_ALL = _${PYTHON_MODN}.so
endif

ifeq ($(USE_MKL),0)
	MKL =
endif

ifeq ($(USE_CUDA),0)
	COMPILE_CUDA =
	CUDA =
	CUDA_WRAP =
	CUDA_OBJECT =
else
	COMPILE_CUDA = cuda
	CUDA_WRAP = ${FORTRAN_CUDA_PATH}/MagTenseCudaBlasICLWrapper.o
	CUDA_OBJECT = ${FORTRAN_CUDA_PATH}/MagTenseCudaBlas.o
endif

ifeq ($(USE_MICROMAG),0)
	MICROMAG =
	MICROMAG_LIB =
else
	MICROMAG = micromagnetism
endif

LIBSRC_SOURCES = ${NUM_INT_PATH}/IOgeneral \
${NUM_INT_PATH}/IntegrationDataTypes \
${NUM_INT_PATH}/quadpack \
${NUM_INT_PATH}/SpecialFunctions \
${TILE_DEMAG_TENSOR_PATH}/TileTensorHelperFunctions \
${TILE_DEMAG_TENSOR_PATH}/TileRectangularPrismTensor \
${TILE_DEMAG_TENSOR_PATH}/TileCircPieceTensor \
${TILE_DEMAG_TENSOR_PATH}/TileCylPieceTensor \
${TILE_DEMAG_TENSOR_PATH}/TilePlanarCoilTensor \
${TILE_DEMAG_TENSOR_PATH}/TileTriangle \
${TILE_DEMAG_TENSOR_PATH}/TileNComponents \
${DEMAG_FIELD_PATH}/DemagFieldGetSolution \
${DEMAG_FIELD_PATH}/spline \
${DEMAG_FIELD_PATH}/MagParameters \
${DEMAG_FIELD_PATH}/IterateMagnetSolution

ifeq ($(USE_MICROMAG),1)
	LIBSRC_SOURCES += ${MICROMAG_PATH}/MicroMagParameters \
	${NUM_INT_PATH}/RKSuitef90 \
	${NUM_INT_PATH}/ODE_Solvers \
	${MICROMAG_PATH}/LLODE_Debug \
	${NUM_INT_PATH}/util \
	${MICROMAG_PATH}/FortranCuda \
	${MICROMAG_PATH}/MagTenseMicroMagPyIO \
	${MICROMAG_PATH}/LandauLifshitzEquationSolver
endif

LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

.PHONY: all clean

all: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${PYTHON_MODN_ALL}

clean:
	cd ${NUM_INT_PATH} && make clean
	cd ${TILE_DEMAG_TENSOR_PATH} && make clean
	cd ${DEMAG_FIELD_PATH} && make clean
	cd ${MICROMAG_PATH} && make clean
	cd ${FORTRAN_CUDA_PATH} && make clean
	${CLEAN_FILES}
	${CLEAN_FOLDER}
	${CLEAN_FOLDER_WIN}

magnetostatic:
	cd ${NUM_INT_PATH} && $(MAKE)
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE)
	cd ${DEMAG_FIELD_PATH} && $(MAKE)

micromagnetism:
	cd ${MICROMAG_PATH} && $(MAKE)

cuda:
	cd ${FORTRAN_CUDA_PATH} && $(MAKE)

${STATIC_LIB}: ${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_LIB} ${DEMAG_FIELD_LIB} ${MICROMAG_LIB} ${CUDA_OBJECT} ${CUDA_WRAP}
	${LIBTOOL}

${PYTHON_MODN_ALL}: ${STATIC_LIB}
	f2py --build-dir . --fcompiler=${FCOMP} ${EXTRA_OPT} -c -m ${PYTHON_MODN} \
	${LIB_OPT} FortranToPythonIO.f90 ${LIBSRC_OBJECTS} $(MKL) $(CUDA)
