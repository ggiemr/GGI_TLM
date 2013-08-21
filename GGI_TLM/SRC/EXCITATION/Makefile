
EXCITATION_OBJS=  \
excitation.o \
excitation_points.o \
excitation_surfaces.o \
mode_excitation.o \
huygens_surface.o 

#set_excitations_in_mesh.o \

default: $(EXCITATION_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o

