
PLOTTING_OBJS=  \
plot_tet_volumes.o \
plot_triangulated_surfaces.o \
plot_line_segments.o \
plot_points.o \
write_volume_list_vtk.o  \
write_surface_list_vtk.o  \
write_line_list_vtk.o  \
write_point_list_vtk.o  \
plot_mesh_volumes.o \
plot_mesh_surfaces.o \
plot_mesh_lines.o \
plot_mesh_points.o \
plot_mesh_boundary.o \
write_volume_mesh_list_vtk.o \
write_surface_mesh_list_vtk.o \
write_line_mesh_list_vtk.o \
write_point_mesh_list_vtk.o \
write_tube_points_vtk.o \
plot_volume_material_cells.o \
plot_surface_material_faces.o

default: $(PLOTTING_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
