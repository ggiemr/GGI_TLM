
MESH_OBJS=  \
reset_mesh_parameters.o \
set_mesh_parameters.o \
general_mesh_subroutines.o \
build_volume_mesh.o \
build_surface_mesh.o \
build_line_mesh.o \
build_point_mesh.o \
mesh_triangle.o \
mesh_tet.o \
mesh_line_segment.o \
read_mesh.o \
write_mesh.o \
trim_mesh.o \
mesh_triangle_new.o \
mesh_tet_new.o \
local_grid_subroutines.o

#model_builder_subroutines.o \

default: $(MESH_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
