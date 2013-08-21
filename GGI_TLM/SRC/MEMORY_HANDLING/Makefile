
MEMORY_HANDLING_OBJS=  \
memory_stuff.o 

default: $(MEMORY_HANDLING_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
