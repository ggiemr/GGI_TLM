PCB_SIMULATION_SETUP_OBJS=  \
process_gerber.o

default: $(PCB_SIMULATION_SETUP_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
