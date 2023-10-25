#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran
USE_CVODE =
CVODE_ROOT =
USE_MATLAB =
MATLAB_INCLUDE =
FFLAGS =
#=======================================================================
#       List all source files required for the project
#=======================================================================
LIBSRC_SOURCES = IOgeneral IntegrationDataTypes quadpack \
SpecialFunctions NumericalIntegrals RKSuitef90 ODE_Solvers util

LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

ifeq ($(OS),Windows_NT)
	LIBTOOL = LIB $? /OUT:$@
	LIB_SUFFIX = .lib
	CLEAN_FILES = -del *.lib *.o *.mod
	
else
	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	LIB_SUFFIX = .a
	CLEAN_FILES = -rm -f .a *.o *.mod
endif
STATIC_LIB = libNumericalIntegration${LIB_SUFFIX}

ifeq ($(USE_CVODE),0)
	CVODE_OPT =
else
	CVODE_OPT = -I"${CVODE_ROOT}"/fortran
endif

ifeq ($(USE_MATLAB),0)
	MATLAB_OPT =
else
	MATLAB_OPT = -I"${MATLAB_INCLUDE}"
endif
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .f90
.PHONY: all clean

all: ${STATIC_LIB}

clean:
	${CLEAN_FILES}

.f90.o:
	$(FC) $(FFLAGS) ${CVODE_OPT} ${MATLAB_OPT} -c $< -o $@

${STATIC_LIB}: ${LIBSRC_OBJECTS}
	${LIBTOOL}
