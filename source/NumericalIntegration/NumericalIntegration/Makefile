#=======================================================================
#                       compiler names and flags
#=======================================================================
# FC = ifort
FC = gfortran

#=======================================================================
#       List all source files required for the project
#=======================================================================
# names (without suffix), f90 sources
LIBSRC_SOURCES = IOgeneral IntegrationDataTypes quadpack \
SpecialFunctions NumericalIntegrals RKSuitef90 ODE_Solvers util

LIB_NAME = libNumericalIntegration

ifeq ($(OS),Windows_NT)
	ifeq (${FC}, ifort)
		FFLAGS = /nologo /real-size:64 /O3 /assume:nocc_omp /Qopenmp /fpp \
		/fpe:0 /fp:source /fp:precise /libs:static /DUSE_CVODE=0
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 \
		-cpp -DUSE_CVODE=0
	endif

	LIBTOOL = LIB $? /OUT:$@
	STATIC_LIB = ${LIB_NAME}.lib
	CLEAN_FILES = -del ${STATIC_LIB} *.o *.mod
	
else
	ifeq (${FC}, ifort)
		FFLAGS = -r8 -O3 -assume nocc_omp -qopenmp -fpp -fpe0 -fp-model source \
		-fp-model precise -fpic -libs:static -DUSE_CVODE=0
	else ifeq (${FC}, gfortran)
		FFLAGS = -fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512 \
		-cpp -DUSE_CVODE=0
	endif

	ifeq (${UNAME}, Darwin)
  		LIBTOOL = libtool -static -o $@ $?
	else
  		LIBTOOL = ar src $@ $?
	endif
	STATIC_LIB = ${LIB_NAME}.a
	CLEAN_FILES = -rm -f ${STATIC_LIB} *.o *.mod
endif

#=======================================================================
#						Relevant suffixes
#=======================================================================
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

.SUFFIXES: .for .f90

#=======================================================================
#							Targets
#=======================================================================
.PHONY: all clean

all: ${STATIC_LIB}

clean:
	${CLEAN_FILES}

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@
.for.o:
	$(FC) $(FFLAGS) -c $< -o $@

${STATIC_LIB}: ${LIBSRC_OBJECTS}
	${LIBTOOL}
