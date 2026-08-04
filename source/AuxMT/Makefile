#=======================================================================
#                       compiler names and flags
#=======================================================================
FC = ifx
FFLAGS =
USE_CVODE =
CVODE_ROOT =
USE_MATLAB =
MATLAB_INCLUDE =
#=======================================================================
#       List all source files required for the project
#=======================================================================

LIBSRC_SOURCES = sort_mod omp_mod omp_timer_mod omp_lock_mod timer_mod trace_mod

LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	OPT = -I${CONDA_PREFIX}/opt/compiler/include/intel64
	LIBTOOL = LIB $? /OUT:$@
	LIB_SUFFIX = .lib
	RM = Remove-Item -Force -ErrorAction SilentlyContinue
	
else
	OPT =
	LIB_SUFFIX = .a
	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
endif

STATIC_LIB = libAuxMT${LIB_SUFFIX}

ifeq ($(USE_CVODE),0)
	CVODE_OPT =
else
	CVODE_OPT = -I${CVODE_ROOT}/fortran
endif

ifeq ($(USE_MATLAB),0)
	MATLAB_OPT =
else
	MATLAB_OPT = -I${MATLAB_INCLUDE}
endif

#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: ${STATIC_LIB}

clean:
	-$(RM) *${LIB_SUFFIX}
	-$(RM) *.o
	-$(RM) *.mod

.f90.o:
	$(FC) $(FFLAGS) ${OPT} ${CVODE_OPT} ${MATLAB_OPT} -c $< -o $@

${STATIC_LIB}: ${LIBSRC_OBJECTS}
	${LIBTOOL}
