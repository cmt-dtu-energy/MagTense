#==========================
# Compiler names and flags
#==========================
USE_CUDA = 1
USE_CVODE = 0
USE_MICROMAG = 1
CPP = icx

ifeq ($(OS),Windows_NT)
	FC = ifort
else
	ifeq (${UNAME}, Darwin)
  		FC = gfortran
	else
  		FC = ifort
	endif
endif

ifeq (${FC}, gfortran)
	USE_CUDA = 0
	USE_MKL = 0
	USE_MICROMAG = 0
endif

ifeq (${FC}, ifort)
	ifeq ($(OS),Windows_NT)
		FFLAGS = /nologo /real-size:64 /O3 /assume:nocc_omp /Qopenmp /fpp \
		/fpe:0 /fp:source /fp:precise /libs:static /DUSE_CVODE=${USE_CVODE} /DUSE_MATLAB=0 \
		/DUSE_CUDA=${USE_CUDA} /DUSE_MICROMAG=${USE_MICROMAG}
	else
		FFLAGS = -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source \
		-fp-model precise -fpic -libs:static -DUSE_CVODE=${USE_CVODE} -DUSE_MATLAB=0 \
		-DUSE_CUDA=${USE_CUDA} -DUSE_MICROMAG=${USE_MICROMAG}
	endif
else ifeq (${FC}, gfortran)
FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 -cpp \
-DUSE_MICROMAG=0 -DUSE_CVODE=${USE_CVODE}
endif

PYTHON_MODN = magtensesource

NUM_INT_PATH = ../../../../source/NumericalIntegration/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = ../../../../source/TileDemagTensor/TileDemagTensor
DEMAG_FIELD_PATH = ../../../../source/DemagField/DemagField
MICROMAG_PATH = ../../../../source/MagTenseMicroMag
FORTRAN_CUDA_PATH = ../../../../source/MagTenseFortranCuda/cuda

VPATH = ${NUM_INT_PATH}:${TILE_DEMAG_TENSOR_PATH}:${DEMAG_FIELD_PATH}:\
${MICROMAG_PATH}:${FORTRAN_CUDA_PATH}

ifeq ($(OS),Windows_NT)
	MKL_ROOT = ${CONDA_PREFIX}/Library
	MKL_ONEAPI = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
	CUDA = -L${CONDA_PREFIX}/Lib/x64 -lcublas -lcudart -lcuda -lcusparse
	MKL = -L${MKL_ONEAPI}/lib/intel64 -lmkl_rt -lmkl_blas95_lp64
	MKL_OMP = -L${MKL_ONEAPI}/windows/compiler/lib/intel64_win -llibiomp5md
	CVODE_ROOT = "C:\Program Files (x86)\sundials-4.1.0\instdir"

	ifeq (${FC}, ifort)
		FCOMP = intelvem
		EXTRA_OPT = -DNO_APPEND_FORTRAN -DUPPERCASE_FORTRAN --no-lower \
		--opt='${FFLAGS} /assume:nounderscore /names:uppercase'
		OPT = ${MKL_ROOT}/include -I${MKL_ONEAPI}/include \
		-I${MKL_ROOT}/include/intel64/lp64 -I${CVODE_ROOT}/fortran \
		-I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} -I${DEMAG_FIELD_PATH} \
		-I${MICROMAG_PATH} -I${FORTRAN_CUDA_PATH}
		LIB_OPT = -I:${OPT} -L. -llibsrc
	else ifeq (${FC}, gfortran)
		FCOMP = gnu95 --compiler=mingw32
		EXTRA_OPT = --opt='${FFLAGS}'
		LIB_OPT = -L. -llibsrc
	endif
		
	LIBTOOL = LIB $? /OUT:$@
	LIB_SUFFIX = .lib
	CLEAN_FILES = -del *.lib *.o *.mod *.pyd
	CLEAN_FOLDER = -rmdir /s /q Release 2>NUL || VER>NUL
	CLEAN_FOLDER_WIN = -for /f %%i in ('dir /a:d /b src.win-*') do rd /s /q %%i
	PYTHON_MODN_ALL = _${PYTHON_MODN}.dll

else
	MKL_ROOT = ${CONDA_PREFIX}
	CVODE_ROOT = /usr/local/sundials-4.1.0/instdir

	ifeq (${FC}, ifort)
		FCOMP = intelem
		OPT = ${MKL_ROOT}/include -I${MKL_ROOT}/include/intel64/lp64 \
		-I${CVODE_ROOT}/fortran -I${NUM_INT_PATH} -I${TILE_DEMAG_TENSOR_PATH} \
		-I${DEMAG_FIELD_PATH} -I${MICROMAG_PATH} -I${FORTRAN_CUDA_PATH}
		LIB_OPT = -I${OPT} -L. -lsrc
	else ifeq (${FC}, gfortran)
		FCOMP = gnu95
		LIB_OPT = -L. -lsrc
	endif

	EXTRA_OPT = --opt='${FFLAGS}'
	CUDA = -L${CONDA_PREFIX}/lib -lcublas -lcudart -lcuda -lcusparse
	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -lpthread -lm -ldl -liomp5 -lmkl_blas95_lp64

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	
	LIB_SUFFIX = .a
	CLEAN_FILES = -rm *.a *.o *.mod *.so 
	CLEAN_FOLDER = -rm -r src.*
	CLEAN_FOLDER_WIN = -true
	PYTHON_MODN_ALL = _${PYTHON_MODN}.so
endif

STATIC_LIB = libsrc${LIB_SUFFIX}
NUM_INT_LIB = ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX}
TILE_DEMAG_TENSOR_LIB = ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX}
DEMAG_FIELD_LIB = ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX}
MICROMAG_LIB = ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX}

CVODE = -L${CVODE_ROOT}/lib -lsundials_nvecserial -lsundials_sunmatrixdense \
-lsundials_sunlinsoldense -lsundials_fnvecserial_mod -lsundials_cvode \
-lsundials_fsunnonlinsolfixedpoint_mod

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

ifeq ($(USE_CVODE),0)
	CVODE =
	CVODE_OPT =
else
	CVODE_OPT = -I${CVODE_ROOT}/fortran
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
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)' CVODE_OPT=$(CVODE_OPT)
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)'
	cd ${DEMAG_FIELD_PATH} && $(MAKE) FC=$(FC) FFLAGS='$(FFLAGS)' USE_MATLAB=0

micromagnetism:
	cd ${MICROMAG_PATH} && $(MAKE) FFLAGS='$(FFLAGS)' MKL_ROOT=$(MKL_ROOT) USE_MATLAB=0

cuda:
	cd ${FORTRAN_CUDA_PATH} && $(MAKE) CPP=$(CPP)

${STATIC_LIB}: ${NUM_INT_LIB} ${TILE_DEMAG_TENSOR_LIB} ${DEMAG_FIELD_LIB} ${MICROMAG_LIB} ${CUDA_OBJECT} ${CUDA_WRAP}
	${LIBTOOL}

${PYTHON_MODN_ALL}: ${STATIC_LIB}
	f2py --build-dir . --fcompiler=${FCOMP} ${EXTRA_OPT} -c -m ${PYTHON_MODN} \
	${LIB_OPT} FortranToPythonIO.f90 ${LIBSRC_OBJECTS} $(MKL) $(CUDA) $(CVODE)
