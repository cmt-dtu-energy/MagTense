#=======================================================================
#                       compiler names and flags
#=======================================================================
USE_CUDA = 1
USE_CVODE = 0
USE_MICROMAG = 1
USE_MATLAB = 0
MATLAB_INCLUDE ?=
USE_FMM3D ?= 0

# /usr/local/MATLAB/<version>/extern/include (Linux)
# "C:\Program Files\MATLAB\<version>\extern\include" (Win)
CPP = icx
FC = ifx
MKFILE_PATH := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
# The CI job unpacks the sundials artifact into <repo>/cvode, so that is the
# default. A local sundials install lives somewhere else, so allow both the
# environment and the command line to override it, e.g.
#   make ... USE_CVODE=1 CVODE_ROOT="C:/Program Files (x86)/sundials-7.2.1"
CVODE_ROOT ?= ${MKFILE_PATH}/cvode

#=======================================================================
#                    Recursive make on Windows
#
# Recipes are run through conda's msys sh, which mounts ${CONDA_PREFIX}/Library
# as / and binds /bin to Library/usr/bin. The PATH entry Library/bin - the only
# directory holding make.exe - is therefore translated to /bin, where it is
# shadowed, and a recursive make invoked by name fails with
#   /usr/bin/sh: line 1: make: command not found
# Pinning MAKE to the absolute path of the conda make.exe makes every sub-make
# resolve, no matter how the top-level make was invoked.
#=======================================================================
ifeq ($(OS),Windows_NT)
	CONDA_MAKE := $(wildcard $(subst \,/,${CONDA_PREFIX})/Library/bin/make.exe)
	ifneq ($(CONDA_MAKE),)
		MAKE := $(CONDA_MAKE)
	endif

	# The include paths are passed to the compiler unquoted, and they have to
	# stay that way: a quote is what makes make hand the recipe to msys sh
	# instead of exec'ing it, and that shell rewrites every /-leading argument
	# of a native program into a path, so /O3, /fpp and the /D defines silently
	# turn into C:/.../Library/O3 and friends and the build comes out wrong.
	# Unquoted, a path with a space is split into two arguments instead, so the
	# space is taken out of the path rather than quoted around, by asking for
	# the 8.3 short form. CI passes paths that have no spaces and never gets here.
	SPACE :=
	SPACE := $(SPACE) $(SPACE)
	short_path = $(if $(findstring $(SPACE),$1),$(or $(shell cygpath -d -m "$1"),$1),$1)
	ifneq ($(MATLAB_INCLUDE),)
		override MATLAB_INCLUDE := $(call short_path,$(MATLAB_INCLUDE))
	endif
	override CVODE_ROOT := $(call short_path,$(CVODE_ROOT))
endif

#=======================================================================
#                    FMM3D integration (upstream Makefile)
#=======================================================================
FMM3D_DEBUG ?= 0
# Submodule location
FMM3D_DIR      ?= external/FMM3D

# Absolute paths to headers and static archive
FMM3D_ROOT     := $(abspath $(FMM3D_DIR))
FMM3D_LIB      := $(FMM3D_ROOT)/local

.PHONY: fmm3d
fmm3d:
ifeq ($(USE_FMM3D),1)
	@echo "==> FMM3D: building via upstream makefile (install)"
	@cd "$(FMM3D_DIR)" && \
	  $(MAKE) install PREFIX=$(abspath $(FMM3D_DIR)/local) DO_DEBUG=$(FMM3D_DEBUG) FAST_KER=OFF
else
	@echo "USE_FMM3D=0 -> skipping FMM3D build"
endif

# Make high-level targets depend on FMM only when enabled
ifeq ($(USE_FMM3D),1)
ALL_DEPS := fmm3d
PY_DEPS  := fmm3d
else
ALL_DEPS :=
PY_DEPS  :=
endif
#=======================================================================

ifeq (${UNAME}, Darwin)
	FC = gfortran
	USE_CUDA = 0
	USE_MKL = 0
	USE_MICROMAG = 0
	USE_MATLAB = 0
endif

ifeq (${FC}, ifx)
	ifeq ($(OS),Windows_NT)
		# /heap-arrays and /traceback match the Linux flags below. Without /heap-arrays
		# ifx puts array temporaries on the stack, and create_CSR_matrix builds three of
		# them at once - pack(rows), pack(columns), pack(values) - sized by the length of
		# the interpolation stencil list. That list grows as ~51x the number of tiles, so
		# from roughly 15000 tiles the three temporaries exceed the stack and the process
		# dies with STATUS_STACK_OVERFLOW (0xC00000FD) inside computeDifferentialOperators-
		# FromMesh_DirectLap. /traceback is what makes any such abort say where it happened.
		FFLAGS = /O3 /fpp /real-size:64 /Qopenmp /assume:nocc_omp /fpe:0 \
			/heap-arrays:1024 /traceback \
			/fp:source /nologo /DUSE_CVODE=${USE_CVODE} /DUSE_MATLAB=${USE_MATLAB} \
			/DUSE_CUDA=${USE_CUDA} /DUSE_MICROMAG=${USE_MICROMAG} /DUSE_FMM3D=${USE_FMM3D}
	else
		FFLAGS = -O3 -fpp -real-size 64 -qopenmp -assume nocc_omp -fpe0 \
			-heap-arrays 1024 -traceback \
			-fp-model=source -fpic -nologo -DUSE_CVODE=${USE_CVODE} \
			-DUSE_MATLAB=${USE_MATLAB} -DUSE_CUDA=${USE_CUDA} \
			-DUSE_MICROMAG=${USE_MICROMAG} -DUSE_FMM3D=${USE_FMM3D}

	endif
else ifeq (${FC}, gfortran)
	FFLAGS = -O3 -fdefault-real-8 -fopenmp -ffree-line-length-512 -cpp -fPIC \
		-DUSE_MICROMAG=0 -DUSE_CVODE=${USE_CVODE} -DUSE_FMM3D=${USE_FMM3D}

endif

PYTHON_MODN = magtensesource
PYTHON_LIBPATH = python/src/magtense/lib

AUXMT_PATH = source/AuxMT
NUM_INT_PATH = source/NumericalIntegration
TILE_DEMAG_TENSOR_PATH = source/TileDemagTensor
DEMAG_FIELD_PATH = source/DemagField
MICROMAG_PATH = source/MagTenseMicroMag
FORTRAN_CUDA_PATH = source/MagTenseFortranCuda/cuda
STANDALONE_PATH = source/MagTense_StandAlone
FORCEINTEGRATOR_PATH = source/MagneticForceIntegrator

VPATH = ${AUXMT_PATH}:${NUM_INT_PATH}:${TILE_DEMAG_TENSOR_PATH}:${DEMAG_FIELD_PATH}:\
${MICROMAG_PATH}:${FORTRAN_CUDA_PATH}:${STANDALONE_PATH}:${FORCEINTEGRATOR_PATH}

ifeq ($(OS),Windows_NT)
	CONDA_PATH = $(subst \,/,${CONDA_PREFIX})
	CUDA_ROOT = ${CONDA_PATH}/Library/lib
	MKL = -L${CONDA_PATH}/Library/lib -lmkl_intel_lp64_dll -lmkl_intel_thread_dll \
		-lmkl_core_dll -lmkl_blas95_lp64 -llibiomp5md
		
	AUXMT_PATH_ROOT     := $(abspath $(AUXMT_PATH))
	
	ifeq ($(USE_FMM3D),0)	
		LDFLAGS = '/DEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcmt.lib /LIBPATH:${CONDA_PATH}/Library/lib'
	else
		LDFLAGS = '/DEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcmt.lib /LIBPATH:${CONDA_PATH}/Library/lib /LIBPATH:${FMM3D_LIB}'
	endif
	
	LIB_SUFFIX = .lib
	PY_MOD_SUFFIX = .pyd
	CVODE_SUFFIX = _static

	ifeq (${FC}, ifx)
		EXTRA_FFLAGS = "${FFLAGS} /assume:underscore /names:lowercase"
		OPT = ${CONDA_PATH}/Library/include \
			-I${CONDA_PATH}/Library/include/intel64/lp64 \
			-I${CONDA_PATH}/opt/compiler/include/intel64
		LIB_OPT = -llibDemagField -llibTileDemagTensor -llibNumericalIntegration
	endif
else
 	MKL = -L${CONDA_PREFIX}/lib -lmkl_rt -liomp5 -lmkl_blas95_lp64 -lpthread -lm -ldl
	CUDA_ROOT = ${CONDA_PREFIX}/lib
# adding LD flags here to aviod "cannot enable executable stack as shared object requires" on newer Linux distributions.
	LDFLAGS = "-Wl,-z,noexecstack"
#	LDFLAGS += '-lstdc++ -liomp5'
	LIB_SUFFIX = .a
	PY_MOD_SUFFIX = .so
	CVODE_SUFFIX =

	ifeq (${FC}, ifx)
		EXTRA_FFLAGS = "${FFLAGS}"
		OPT = ${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/intel64/lp64
		LIB_OPT = -lNumericalIntegration -lTileDemagTensor -lDemagField
	endif
endif

ifeq (${FC}, gfortran)
	LIB_OPT =
endif

ifeq ($(USE_MKL),0)
	MKL =
endif


AUXMT = auxmt
ifeq ($(OS),Windows_NT)
	LIB_OPT += -llibAuxMT
else
	LIB_OPT += -lAuxMT
endif

ifeq ($(USE_MICROMAG),0)
	MICROMAG =
	ifeq ($(OS),Windows_NT)
		CP_LIB = cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
	else
		CP_LIB = cp ${AUXMT_PATH}/libAuxMT${LIB_SUFFIX} .
		CP_LIB += && cp ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX} .
		CP_LIB += && cp ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX} .
		CP_LIB += && cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
	endif
else
	MICROMAG = micromagnetism
	ifeq ($(OS),Windows_NT)
		LIB_OPT = -llibAuxMT
		CP_LIB = cp ${AUXMT_PATH}/libAuxMT${LIB_SUFFIX} .
		
		LIB_OPT += -llibMagTenseMicroMag
		#CP_LIB = cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
		CP_LIB += && cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
		
	else
		LIB_OPT += -lMagTenseMicroMag
		CP_LIB = cp ${AUXMT_PATH}/libAuxMT${LIB_SUFFIX} .
		CP_LIB += && cp ${NUM_INT_PATH}/libNumericalIntegration${LIB_SUFFIX} .
		CP_LIB += && cp ${TILE_DEMAG_TENSOR_PATH}/libTileDemagTensor${LIB_SUFFIX} .
		CP_LIB += && cp ${DEMAG_FIELD_PATH}/libDemagField${LIB_SUFFIX} .
		CP_LIB += && cp ${MICROMAG_PATH}/libMagTenseMicroMag${LIB_SUFFIX} .
	endif



endif
#===================== FMM3D Integration ===========================
ifeq ($(USE_FMM3D),0)
  FMM3D =
else
  ifeq ($(OS),Windows_NT)
    # 1. Use the standard linker flag (-lfmm3d looks for fmm3d.lib)
    FMM3D = -L"${FMM3D_LIB}" -lfmm3d
    
    # 2. Use 'cp' (since your shell is sh) and RENAME the file during copy
    # This creates fmm3d.lib in your root so the linker finds it easily
    CP_LIB += && cp "${FMM3D_LIB}/libfmm3d.lib" ./fmm3d.lib
  else
    FMM3D = -L${FMM3D_LIB} -lfmm3d
    LDFLAGS = "-Wl,-z,noexecstack -Wl,-rpath,${FMM3D_LIB}"
    CP_LIB += && cp ${FMM3D_LIB}/libfmm3d${LIB_SUFFIX} .
  endif
endif
#===================================================================



ifeq ($(USE_MATLAB),0)
	FORCEINTEGRATOR =
else
	FORCEINTEGRATOR = forceintegrator
	CP_LIB += && cp ${FORCEINTEGRATOR_PATH}/libMagneticForceIntegrator${LIB_SUFFIX} .
endif

ifeq ($(USE_CUDA),0)
	COMPILE_CUDA =
	CUDA =
else
	COMPILE_CUDA = cuda
	CUDA = -L${CUDA_ROOT} -lcublas -lcudart -lcusparse
	CP_LIB += && cp ${FORTRAN_CUDA_PATH}/libCuda${LIB_SUFFIX} .
	ifeq ($(OS),Windows_NT)
		CUDA += -lcuda
		LIB_OPT += -llibCuda
	else
		LIB_OPT += -lCuda
	endif
endif

ifeq ($(USE_CVODE),0)
	CVODE =
else
	CVODE = -L${CVODE_ROOT}/lib -lsundials_core${CVODE_SUFFIX} -lsundials_cvode${CVODE_SUFFIX} \
		-lsundials_fcore_mod -lsundials_fcvode_mod${CVODE_SUFFIX} -lsundials_fnvecserial_mod${CVODE_SUFFIX} \
		-lsundials_fsunmatrixdense_mod${CVODE_SUFFIX} -lsundials_fsunlinsolspgmr_mod${CVODE_SUFFIX}
endif

INCLUDE_OBJ = ${MKFILE_PATH}/${AUXMT_PATH} \
	-I${MKFILE_PATH}/${NUM_INT_PATH} \
	-I${MKFILE_PATH}/${TILE_DEMAG_TENSOR_PATH} \
	-I${MKFILE_PATH}/${DEMAG_FIELD_PATH} \
	-I${MKFILE_PATH}/${MICROMAG_PATH} \
	-I${MKFILE_PATH}/${FORTRAN_CUDA_PATH} 

PYTHON_MODN_ALL = _${PYTHON_MODN}${PY_MOD_SUFFIX}

#=======================================================================
#                    Build configuration consistency
#
# The USE_* settings are baked into the objects as preprocessor defines, so the
# static libraries and the python extension have to be built with the same set.
# Make cannot tell that a flag changed - the sources are unmodified, so nothing
# is recompiled - and a mismatch therefore surfaces only as a wall of unresolved
# symbols when the extension is linked. Building the libraries with USE_FMM3D=1
# and then linking with USE_FMM3D=0 leaves every FMM3D symbol undefined, and the
# reverse leaves the FMM3D entry points out of the extension.
#
# The settings are therefore recorded when the libraries are built and compared
# before the extension is linked, so that a mismatch is reported in terms of the
# flags that caused it.
#=======================================================================
BUILD_FLAGS_FILE := .build_flags
BUILD_FLAGS := USE_CUDA=${USE_CUDA} USE_CVODE=${USE_CVODE} USE_MATLAB=${USE_MATLAB} USE_MICROMAG=${USE_MICROMAG} USE_FMM3D=${USE_FMM3D}

# Written by the library targets once they have succeeded
RECORD_FLAGS = @echo "${BUILD_FLAGS}" > ${BUILD_FLAGS_FILE}

.PHONY: check-flags
check-flags:
	@if [ -f ${BUILD_FLAGS_FILE} ] && [ "`cat ${BUILD_FLAGS_FILE}`" != "${BUILD_FLAGS}" ]; then \
		echo "ERROR: the libraries and this link step were configured differently."; \
		echo "       libraries built with: `cat ${BUILD_FLAGS_FILE}`"; \
		echo "       linking with:         ${BUILD_FLAGS}"; \
		echo "       Re-run both steps with the same USE_* settings, or 'make clean' first."; \
		exit 1; \
	fi

#=======================================================================
#                    Optional dependency locations
#
# USE_CVODE=1 needs a sundials install and USE_MATLAB=1 needs MATLAB's extern
# headers, neither of which the repository carries. CI supplies both - it
# unpacks the sundials artifact into <repo>/cvode and passes MATLAB_INCLUDE on
# the command line - so a local build that leaves them unset is the only way to
# reach these paths, and the failures are unrecognisable without this check:
# a missing CVODE_ROOT surfaces as "error #7002" on every sundials module, and
# an empty MATLAB_INCLUDE expands to a bare -I that swallows the following -c,
# turning the compile into a link and burying the build under
# "libifcoremt.lib(for_main.obj) : error LNK2019: unresolved external MAIN__".
#=======================================================================
.PHONY: check-config
check-config:
ifeq ($(OS),Windows_NT)
	@case "${CVODE_ROOT}${MATLAB_INCLUDE}" in *" "*) \
		echo "ERROR: a dependency path still contains a space after 8.3 shortening."; \
		echo "       CVODE_ROOT     = ${CVODE_ROOT}"; \
		echo "       MATLAB_INCLUDE = ${MATLAB_INCLUDE}"; \
		echo "       The compiler is invoked without quotes - see the note at the top"; \
		echo "       of this Makefile - so the path is split on the space. Either move"; \
		echo "       the dependency somewhere without spaces, or enable 8.3 names on"; \
		echo "       the volume (fsutil 8dot3name set 0) and recreate the directory."; \
		exit 1;; \
	esac
endif
ifeq ($(USE_CVODE),1)
	@if [ ! -d "${CVODE_ROOT}/fortran" ]; then \
		echo "ERROR: USE_CVODE=1 but no sundials Fortran modules under CVODE_ROOT."; \
		echo "       CVODE_ROOT = ${CVODE_ROOT}"; \
		echo "       Point it at your sundials install, e.g."; \
		echo "         make ... USE_CVODE=1 CVODE_ROOT=\"C:/Program Files (x86)/sundials-7.2.1\""; \
		exit 1; \
	fi
endif
ifeq ($(USE_MATLAB),1)
	@if [ -z "${MATLAB_INCLUDE}" ]; then \
		echo "ERROR: USE_MATLAB=1 requires MATLAB_INCLUDE to be set."; \
		echo "       Point it at MATLAB's extern/include, e.g."; \
		echo "         make ... USE_MATLAB=1 MATLAB_INCLUDE=\"C:/Program Files/MATLAB/R2024b/extern/include\""; \
		exit 1; \
	fi
	@if [ ! -f "${MATLAB_INCLUDE}/fintrf.h" ]; then \
		echo "ERROR: no fintrf.h under MATLAB_INCLUDE."; \
		echo "       MATLAB_INCLUDE = ${MATLAB_INCLUDE}"; \
		exit 1; \
	fi
endif

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: $(ALL_DEPS) ${AUXMT} magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} 

standalone: magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${FORCEINTEGRATOR} standalone

python: $(PY_DEPS) ${AUXMT} magnetostatic ${MICROMAG} ${COMPILE_CUDA} ${PYTHON_MODN_ALL}

python-win: $(PY_DEPS) ${PYTHON_MODN_ALL}

clean:
	cd ${AUXMT_PATH} && ${MAKE} clean
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	cd $(STANDALONE_PATH) && ${MAKE} clean
	cd $(FORCEINTEGRATOR_PATH) && ${MAKE} clean
	rm -f *${LIB_SUFFIX} *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/*${LIB_SUFFIX} ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build
	rm -f ${BUILD_FLAGS_FILE}

clean_full:
	cd ${FMM3D_DIR} && ${MAKE} clean
	cd ${AUXMT_PATH} && ${MAKE} clean
	cd ${NUM_INT_PATH} && ${MAKE} clean
	cd ${TILE_DEMAG_TENSOR_PATH} && ${MAKE} clean
	cd ${DEMAG_FIELD_PATH} && ${MAKE} clean
	cd ${MICROMAG_PATH} && ${MAKE} clean
	cd ${FORTRAN_CUDA_PATH} && ${MAKE} clean
	cd $(STANDALONE_PATH) && ${MAKE} clean
	cd $(FORCEINTEGRATOR_PATH) && ${MAKE} clean
	rm -f *${LIB_SUFFIX} *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/*${LIB_SUFFIX} ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build
	rm -f ${BUILD_FLAGS_FILE}

auxmt: check-config
	cd ${AUXMT_PATH} && ${MAKE} FC=${FC} FFLAGS="${FFLAGS}" USE_CVODE=${USE_CVODE} CVODE_ROOT="${CVODE_ROOT}" USE_MATLAB=${USE_MATLAB} MATLAB_INCLUDE="${MATLAB_INCLUDE}"
	${RECORD_FLAGS}

clean-build:
	rm -f ${PYTHON_LIBPATH}/*${PY_MOD_SUFFIX}
	rm -rf ${PYTHON_LIBPATH}/build

magnetostatic: check-config
	cd ${NUM_INT_PATH} && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT="$(CVODE_ROOT)" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="$(MATLAB_INCLUDE)"
	cd ${TILE_DEMAG_TENSOR_PATH} && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT="$(CVODE_ROOT)" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="$(MATLAB_INCLUDE)"
	cd ${DEMAG_FIELD_PATH}  && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT="$(CVODE_ROOT)" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="$(MATLAB_INCLUDE)"
	${RECORD_FLAGS}

micromagnetism: check-config
	cd ${MICROMAG_PATH} && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' USE_CVODE=$(USE_CVODE) CVODE_ROOT="$(CVODE_ROOT)" USE_MATLAB=$(USE_MATLAB) MATLAB_INCLUDE="$(MATLAB_INCLUDE)"
	${RECORD_FLAGS}

cuda:
	cd ${FORTRAN_CUDA_PATH} && $(MAKE) CPP=$(CPP)

forceintegrator: check-config
	cd $(FORCEINTEGRATOR_PATH) && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}' MATLAB_INCLUDE="$(MATLAB_INCLUDE)"

standalone:
	cd $(STANDALONE_PATH) && $(MAKE) FC=$(FC) FFLAGS='${FFLAGS}'
	mkdir build
	cp $(STANDALONE_PATH)/MagTense.x build/MagTense.x

info:
	@echo Using Fortran compiler: $(FC)
	@echo Fortran flags: $(FFLAGS)
	@echo Linker flags: $(LDFLAGS)
	@echo MKL flags: $(MKL)
	@echo CUDA enabled: $(USE_CUDA)
	@echo CVODE enabled: $(USE_CVODE)
	@echo Micromagnetics enabled: $(USE_MICROMAG)
	@echo FMM3D enabled: $(USE_FMM3D)
	@echo MATLAB enabled: $(USE_MATLAB)
	@echo MATLAB include: $(MATLAB_INCLUDE)
	@echo CVODE root: $(CVODE_ROOT)
	@echo Include paths: -I$(INCLUDE_OBJ)
	@echo Libraries: -L${MKFILE_PATH} ${LIB_OPT}



.PHONY: test
test:
	@command -v fpm >/dev/null 2>&1 || (echo "ERROR: fpm not found. Install with: conda install -c conda-forge fpm" && exit 1)
	@echo "==> Cleaning all fpm test build directories"
	@find tests -mindepth 2 -maxdepth 2 -type d -name build -exec rm -rf {} +
	@echo "==> Running all fpm test harnesses"
	cd tests && FPM_FC=$(FC) FPM_FFLAGS='$(FFLAGS)' bash run_all_fpm_tests.sh




${PYTHON_MODN_ALL}: check-config check-flags
	${CP_LIB}
	FC=${FC} FFLAGS=${EXTRA_FFLAGS} LDFLAGS=${LDFLAGS} \
		python -m numpy.f2py -c -m ${PYTHON_MODN} \
		--build-dir ${PYTHON_LIBPATH}/build -I${OPT} -I${INCLUDE_OBJ} \
		-L${MKFILE_PATH} ${LIB_OPT} python/FortranToPythonIO.f90 ${MKL} ${CUDA} ${CVODE} ${FMM3D}
	cp *${PY_MOD_SUFFIX} ${PYTHON_LIBPATH}/
