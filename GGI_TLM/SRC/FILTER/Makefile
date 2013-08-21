
FILTER_OBJS=  \
filter_types.o \
filter_functions.o \
filter_operators.o \
filter_subroutines.o 

#
default: $(FILTER_OBJS)
#
%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o


