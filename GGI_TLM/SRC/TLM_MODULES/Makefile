
MODULE_OBJS=  \
Geometry_Modules.o \
TLM_Modules.o \
TLM_material_modules.o \
TLM_excitation_modules.o \
TLM_output_modules.o \
TLM_mode_stir_modules.o \
Cable_Modules.o \
Mesh_Modules.o \
pul_LC_Modules.o \
Filter_fit_modules.o \
File_information.o \
post_process_modules.o 

default: $(MODULE_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
