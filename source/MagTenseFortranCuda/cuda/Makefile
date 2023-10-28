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
CUDA_FILE = $(addsuffix .cu, ${CUDA_SOURCE})
CUDA_WRAP_FILE = $(addsuffix .cxx, ${CUDA_WRAP_SOURCE})
#=======================================================================
#
#=======================================================================
ifeq ($(OS),Windows_NT)
	CUDA_OBJECT = $(addsuffix .obj, ${CUDA_SOURCE})
	CUDA_WRAP = $(addsuffix .obj, ${CUDA_WRAP_SOURCE})
	CPPFLAGS = /Fo:${CUDA_WRAP}
	CLEAN_FILES = -del *.obj
else
	CUDA_OBJECT = $(addsuffix .o, ${CUDA_SOURCE})
	CUDA_WRAP = $(addsuffix .o, ${CUDA_WRAP_SOURCE})
	CUFLAGS = -shared -Xcompiler -fPIC -I ${CONDA_PREFIX}/include
	CPPFLAGS = ${CFLAGS} -o ${CUDA_WRAP}
	CLEAN_FILES = -rm *.o
endif

.PHONY: all clean

all: cuda cuda_wrap

clean:
	${CLEAN_FILES}

cuda: 
	nvcc ${CUFLAGS} -c ${CUDA_FILE} -o ${CUDA_OBJECT}

cuda_wrap: 
	${CPP} ${CPPFLAGS} -c ${CUDA_WRAP_FILE} -o ${CUDA_WRAP}
