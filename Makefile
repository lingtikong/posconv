.SUFFIXES : .o .c .f .f90
#
# Machine dependent info
FC = /opt/intel/bin/ifort
#FC = gfortran
CC = cc

# optimization flags
OFLAG = -O -g -traceback

# Fortran and C flags
FFLAGS = $(DFLAGS) $(OFLAG) $(DEBUG)
CFLAGS = $(DFLAGS) $(OFLAG) $(DEBUG)
 
# Fortran 90/95 format
FREE = -free
# FREE = -ffree-form -ffree-line-length-none
#
# 
BASE = posconv
MAIN = ${BASE}
#
# source and rules
#====================================================================
SOURCE = prec.o vardef.o date.o lattconv.o readpos.o writepos.o \
         posvasp.o pospwscf.o posxyz.o posgrof.o posbgf.o       \
         lmp_atom.o siesta.o abinit.o latgen.o posart.o         \
         posvel.o posms.o lammps.o screen.o operation.o         \
         realign.o nominal.o neighbor.o error.o main.o

all:  ${MAIN}

${MAIN}:  $(SOURCE)
	$(FC) $(FREE) $(SOURCE) $(OFLAG) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${MAIN}

tar:
	rm -f ${BASE}.tar; tar -czvf ${BASE}.tar.gz *.f90 Makefile


.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
