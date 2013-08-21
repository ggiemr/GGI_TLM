
CONSTANTS_OBJS=  \
constants.o 

default: $(CONSTANTS_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
