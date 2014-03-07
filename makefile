#
# default makefile for ifort compiler with more or less
# appropriate options for debugging and high performance
#

# application name
APP = femsim

# list of source files
SRC = globals.f90 read_inputs.f90 model.f90 scattering_factors.f90 fem.f90 femsim.f90 

# list of object files
OBJ = globals.o   read_inputs.o   model.o   scattering_factors.o   fem.o   femsim.o 

# define libraries needed by the linker
#LIBS = -L/export/apps/openmpi_intel_20130618/lib -lmpi
LIBS=
#LIBS = -L/export/apps/openmpi_intel_20130618/lib
# -rpath /state/partition1/apps/intel_20130618/composer_xe_2013.3.163/compiler/lib/

# compiler options for debugging
#FC_DEBUG = ifort -g -debug all -check all -implicitnone -warn all
#FC_DEBUG = ifort -g -debug all -check all -implicitnone -warn all -openmp -fpp
#FC_DEBUG = ifort -g -implicitnone
FC_DEBUG = mpif90 -g -implicitnone -openmp -fpp

# compiler options for optmized running
#FC_OPT = ifort -O3 -xO -ipo -no-prec-div -static
#FC_OPT = mpif90 -O3 -ipo -static
FC_OPT = mpif90 -O2 -openmp -fpp

#CFLAGS=-I/export/apps/openmpi_intel_20130618/include/

# build rules

.SUFFIXES: .f90 .o
.f90.o:
	${FC_OPT} -c $<

debug: ${OBJ}
	${FC_OPT} -o ${APP} ${OBJ} ${LIBS}

opt: ${SRC}
	${FC_OPT} -o ${APP} ${SRC} ${LIBS}

clean:
	rm -f *.mod *.o ${APP}
