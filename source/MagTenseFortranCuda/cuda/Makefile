#==========================
# Compiler names and flags
#==========================
CPP = icx
CUDA_HOST_COMPILER ?=
#=======================================================================
#       List all source files required for the project
#=======================================================================
CUDA_SOURCE = MagTenseCudaBlas
CUDA_EXCH_SOURCE = MagTenseCudaExchange
CUDA_WRAP_SOURCE = MagTenseCudaBlasICLWrapper
#=======================================================================
#                 Relevant suffixes
#=======================================================================
CUDA_FILE = $(CUDA_SOURCE).cu
CUDA_EXCH_FILE = $(CUDA_EXCH_SOURCE).cu
CUDA_WRAP_FILE = $(CUDA_WRAP_SOURCE).cxx

# MagTenseCudaExchange.cu has its own kernel, so it is compiled to machine code for each of these
# architectures (the list dip-fmm uses) plus PTX for the newest. PTX alone is not enough: a driver
# older than the toolkit refuses to JIT it ("PTX was compiled with an unsupported toolchain").
CUDA_ARCHITECTURES ?= 75 80 86 89 90 120
CUDA_GENCODE = $(foreach a,$(CUDA_ARCHITECTURES),-gencode arch=compute_$(a),code=sm_$(a)) \
	-gencode arch=compute_$(lastword $(CUDA_ARCHITECTURES)),code=compute_$(lastword $(CUDA_ARCHITECTURES))

ifeq ($(OS),Windows_NT)
    OBJ_SUFFIX = .obj
    # Use the environment of the surrounding developer shell rather than letting
    # nvcc call vcvars64.bat again, which fails from inside an initialised one.
    NVCC = nvcc
    CUFLAGS = --use-local-env
    CPPFLAGS = /nologo
    CPPFLAGS_OUT = /Fo:
    LIBTOOL = LIB $^ /NOLOGO /OUT:$@
    LIB_SUFFIX = .lib
else
    OBJ_SUFFIX = .o
    # The headers already come from the conda environment, so the compiler has
    # to come from there too. Left to the PATH, a system-wide /usr/local/cuda-*
    # ahead of the environment wins and its nvcc is used against the pinned
    # headers and libraries - an older one then also rejects the environment's
    # gcc outright ("unsupported GNU version"). Fall back to the PATH when the
    # environment has no nvcc of its own.
    NVCC = $(if $(wildcard ${CONDA_PREFIX}/bin/nvcc),${CONDA_PREFIX}/bin/nvcc,nvcc)
    CUFLAGS = -shared -Xcompiler -fPIC -I ${CONDA_PREFIX}/include
    ifneq (${CUDA_HOST_COMPILER},)
        CUFLAGS += -ccbin ${CUDA_HOST_COMPILER}
    endif
    CPPFLAGS = ${CFLAGS}
    CPPFLAGS_OUT = -o 
    LIB_SUFFIX = .a
    ifeq (${UNAME}, Darwin)
        LIBTOOL = libtool -static -o $@ $^
    else
        LIBTOOL = ar src $@ $?
    endif
endif

CUDA_OBJECT = $(CUDA_SOURCE)$(OBJ_SUFFIX)
CUDA_EXCH_OBJECT = $(CUDA_EXCH_SOURCE)$(OBJ_SUFFIX)
CUDA_WRAP   = $(CUDA_WRAP_SOURCE)$(OBJ_SUFFIX)

STATIC_LIB = libCuda$(LIB_SUFFIX)

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean wrap cuda

all: cuda $(STATIC_LIB)

wrap: $(STATIC_LIB)

cuda: $(CUDA_OBJECT) $(CUDA_EXCH_OBJECT)

$(CUDA_OBJECT): $(CUDA_FILE)
	$(NVCC) $(CUFLAGS) -c $(CUDA_FILE) -o $(CUDA_OBJECT)

$(CUDA_EXCH_OBJECT): $(CUDA_EXCH_FILE)
	$(NVCC) $(CUFLAGS) $(CUDA_GENCODE) -c $(CUDA_EXCH_FILE) -o $(CUDA_EXCH_OBJECT)

$(CUDA_WRAP): $(CUDA_WRAP_FILE)
	$(CPP) $(CPPFLAGS) $(CPPFLAGS_OUT)$@ -c $(CUDA_WRAP_FILE)

$(STATIC_LIB): $(CUDA_OBJECT) $(CUDA_EXCH_OBJECT) $(CUDA_WRAP)
	$(LIBTOOL)

clean:
	rm -f *.o *.obj *.a *.lib
