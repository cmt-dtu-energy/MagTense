#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran
CVODE_OPT =
MATLAB_OPT =
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
#=======================================================================
#							Targets
#=======================================================================
.SUFFIXES: .for .f90
.PHONY: all clean

all: ${STATIC_LIB}

clean:
	${CLEAN_FILES}

.f90.o:
	$(FC) $(FFLAGS) ${CVODE_OPT} ${MATLAB_OPT} -c $< -o $@
.for.o:
	$(FC) $(FFLAGS) -c $< -o $@

${STATIC_LIB}: ${LIBSRC_OBJECTS}
	${LIBTOOL}
