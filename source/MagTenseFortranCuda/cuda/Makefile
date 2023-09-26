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

CUDA_OBJECT = $(addsuffix .o, ${CUDA_SOURCE})
CUDA_WRAP = $(addsuffix .o, ${CUDA_WRAP_SOURCE})
#=======================================================================
#
#=======================================================================
ifeq ($(OS),Windows_NT)
	COMPILE_CUDA =
	CPPFLAGS = /Fo:${CUDA_WRAP}
	CLEAN_FILES = -del *.o

else
	COMPILE_CUDA = cuda
	CUFLAGS = -shared -Xcompiler -fPIC -I ${CONDA_PREFIX}/include
	CPPFLAGS = ${CFLAGS} -o ${CUDA_WRAP}
	CLEAN_FILES = -rm *.o
endif

.PHONY: all clean

all: ${COMPILE_CUDA} cuda_wrap

clean:
	${CLEAN_FILES}

cuda: 
	nvcc ${CUFLAGS} -c ${CUDA_FILE} -o ${CUDA_OBJECT}

cuda_wrap: 
	${CPP} ${CPPFLAGS} -c ${CUDA_WRAP_FILE} -o ${CUDA_WRAP}
