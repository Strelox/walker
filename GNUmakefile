# Compiler, compiler options, linker, linker options
FC = gfortran
FCOPTS = -std=f2003 -pedantic -Wall #-fbounds-check
LN = $(FC)
LNOPTS = 
LAPACK =

# Object files
OBJS = accuracy.o io.o config.o random.o randomWalk.o networks.o walker.o


walker: $(OBJS)
	$(LN) $(LNOPTS) -o $@ $^ $(LAPACK)

%.o: %.f90
	$(FC) $(FCOPTS) -c  $<

accuracy.o:

io.o: accuracy.o

config.o: accuracy.o

random.o:

randomWalk.o: accuracy.o io.o random.o

networks.o: accuracy.o io.o random.o

walker.o: accuracy.o io.o random.o


# Clean directory
.PHONY: clean realclean 

clean:
	rm -f *.o *.mod

realclean: clean 
	rm -f walker *.dat
