# Compiler, compiler options, linker, linker options
FC = gfortran
FCOPTS = -std=f2003 -pedantic -Wall -fbounds-check -fall-intrinsics
LN = $(FC)
LNOPTS =

LAPACK95_LIBDIR = /home/frank/lib/LAPACK95
LAPACK95 = lapack95
LAPACK95_MODDIR = /home/frank/lib/LAPACK95/lapack95_modules
LAPACK_LIBDIR = /home/frank/lib/LAPACK
LAPACK = lapack
BLAS_LIBDIR = /home/frank/lib/BLAS
BLAS = blas

# Object files
OBJS = accuracy.o io.o config.o random.o randomWalk.o quantumWalk.o networks.o walker.o


walker: $(OBJS)
	$(LN) $(LNOPTS) -o $@ $^ -L$(LAPACK95_LIBDIR) -L$(LAPACK_LIBDIR) -L$(BLAS_LIBDIR) -l$(LAPACK95) -l$(LAPACK) -l$(BLAS)

%.o: %.f90
	$(FC) $(FCOPTS) -c -I$(LAPACK95_MODDIR) $<

accuracy.o:

io.o: accuracy.o

config.o: accuracy.o

random.o:

randomWalk.o: accuracy.o io.o random.o

quantumWalk.o: accuracy.o io.o random.o

networks.o: accuracy.o io.o random.o

walker.o: accuracy.o io.o random.o


# Clean directory
.PHONY: clean realclean 

clean:
	rm -f *.o *.mod

realclean: clean 
	rm -f walker *.dat
