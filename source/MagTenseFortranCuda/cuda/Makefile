#==========================
# Compiler names and flags
#==========================
CPP = icx
#=======================================================================
#       List all source files required for the project
#=======================================================================
CUDA_SOURCE = MagTenseCudaBlas
CUDA_WRAP_SOURCE = MagTenseCudaBlasICLWrapper
#=======================================================================
#                 Relevant suffixes
#=======================================================================
CUDA_FILE = $(CUDA_SOURCE).cu
CUDA_WRAP_FILE = $(CUDA_WRAP_SOURCE).cxx

ifeq ($(OS),Windows_NT)
    OBJ_SUFFIX = .obj
    CPPFLAGS = /nologo
    CPPFLAGS_OUT = /Fo:
    LIBTOOL = LIB $^ /NOLOGO /OUT:$@
    LIB_SUFFIX = .lib
else
    OBJ_SUFFIX = .o
    CUFLAGS = -shared -Xcompiler -fPIC -I ${CONDA_PREFIX}/include
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
CUDA_WRAP   = $(CUDA_WRAP_SOURCE)$(OBJ_SUFFIX)

STATIC_LIB = libCuda$(LIB_SUFFIX)

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean wrap cuda

all: cuda $(STATIC_LIB)

wrap: $(STATIC_LIB)

cuda:
	nvcc $(CUFLAGS) -c $(CUDA_FILE) -o $(CUDA_OBJECT)

$(CUDA_WRAP):
	$(CPP) $(CPPFLAGS) $(CPPFLAGS_OUT)$@ -c $(CUDA_WRAP_FILE)

$(STATIC_LIB): $(CUDA_OBJECT) $(CUDA_WRAP)
	$(LIBTOOL)

clean:
	rm -f *.o *.obj *.a *.lib