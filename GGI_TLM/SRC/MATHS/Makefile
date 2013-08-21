
MATHS_OBJS=  \
dmatrix.o \
cmatrix.o \
Legendre.o \
bessel.o 

default: $(MATHS_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
