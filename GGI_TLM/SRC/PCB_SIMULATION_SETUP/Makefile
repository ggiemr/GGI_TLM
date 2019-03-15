PCB_SIMULATION_SETUP_OBJS=  \
process_gerber.o  \
add_additional_components.o   \
add_vias.o  \
add_components.o  \
add_dielectric_layers.o  \
circuit_geometry.o  \
write_GGI_TLM_input_file.o  \
write_Spice_input_file.o

default: $(PCB_SIMULATION_SETUP_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
