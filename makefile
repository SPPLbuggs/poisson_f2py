COMPILER = mpif90
COMPFLAG = -fPIC
PETSC_INCL = -I$(PETSC_DIR)/include
PETSC_INCL2 = -I$(PETSC_DIR)/arch-linux2-c-debug/include
LIB_INCL = -L$(PETSC_DIR)/arch-linux2-c-debug/lib
MPICH_LIB = -I/usr/include/mpich
GNU_LIB = -L/usr/lib/x86_64-linux-gnu

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

objects = properties.o laplace_lib.o petsc_lib.o

#--------------------------------------------------------------------------

MPI = $(MPICH_LIB) $(GNU_LIB) -lmpichfort -lmpich
PETSC = $(PETSC_INCL) $(PETSC_INCL2) $(LIB_INCL) -lpetsc

#--------------------------------------------------------------------------
main: $(objects)
	f2py -c -m main main.F90 $(objects) $(MPI) $(PETSC)

%.o: %.F90
	$(COMPILER) $(COMPFLAG) -c $< $(PETSC)
