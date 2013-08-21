
CABLE_MODEL_OBJS=  \
general_cable_subroutines.o  \
build_cell_centre_junction_list.o  \
build_face_junction_list.o  \
create_cable_LCRG_matrices.o  \
create_bundle_LCRG_matrices.o  \
create_bundle_segments.o  \
get_bundle_segment_geometries.o  \
create_cable_meshes.o  \
set_bundle_excitations.o  \
set_bundle_outputs.o  \
write_cable_model_data.o \
read_cable_model_data.o \
create_bundle_cross_section_geometry.o \
check_and_stabilise_cable_model.o \
calculate_TLM_cable_matrices.o \
partition_cable_mesh.o \
cable_model_checks.o

default: $(CABLE_MODEL_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
