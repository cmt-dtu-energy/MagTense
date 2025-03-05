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
	CPPFLAGS = /nologo /Fo:${CUDA_WRAP}
	LIBTOOL = LIB $? /NOLOGO /OUT:$@
	LIB_SUFFIX = .lib
else
	CUFLAGS = -shared -Xcompiler -fPIC -I ${CONDA_PREFIX}/include
	CPPFLAGS = ${CFLAGS} -o ${CUDA_WRAP}
	LIB_SUFFIX = .a
	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
endif

CLEAN_FILES = -rm *.o
STATIC_LIB = libCuda${LIB_SUFFIX}

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: cuda ${STATIC_LIB}

wrap: ${STATIC_LIB}

clean:
	${CLEAN_FILES}

cuda: 
	nvcc ${CUFLAGS} -c ${CUDA_FILE} -o ${CUDA_OBJECT}

${CUDA_WRAP}: 
	${CPP} ${CPPFLAGS} -c ${CUDA_WRAP_FILE} -o $@

${STATIC_LIB}: ${CUDA_OBJECT} ${CUDA_WRAP}
	${LIBTOOL}
