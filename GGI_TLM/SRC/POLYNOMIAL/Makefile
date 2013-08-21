
POLYNOMIAL_OBJS=  \
polynomial_types.o \
polynomial_functions.o \
polynomial_operators.o \
root_finding.o
#
default: $(POLYNOMIAL_OBJS)
#
%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
