
MATERIAL_OBJS=  \
TLM_volume_materials.o \
TLM_surface_materials.o \
output_material_frequency_response.o \
Surface_material_frequency_response.o \
Volume_material_frequency_response.o

default: $(MATERIAL_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o

