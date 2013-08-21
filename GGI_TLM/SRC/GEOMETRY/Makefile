
GEOMETRY_OBJS=  \
geometry_operators.o \
transformations.o \
functions.o \
build_volume_geometry.o \
build_surface_geometry.o \
build_line_geometry.o \
build_point_geometry.o \
get_point.o \
build_surface_sphere.o \
build_surface_rectangular_block.o  \
build_surface_rectangular_block2.o \
build_surface_circle.o \
build_surface_cylinder.o \
build_surface_rectangle.o \
build_surface_xplane.o \
build_surface_yplane.o \
build_surface_zplane.o \
build_surface_quad.o \
build_surface_triangle.o \
build_surface_triangulated_surface.o \
build_surface_vtk_triangulated_surface.o \
build_volume_sphere.o \
build_volume_rectangular_block.o  \
build_volume_rectangular_block2.o \
build_volume_cylinder.o \
build_volume_tet.o \
build_volume_pyramid.o \
build_volume_pyramid_ram.o \
build_volume_tet_mesh.o \
build_line_straight_line.o \
build_line_straight_line2.o \
build_line_arc.o 

default: $(GEOMETRY_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
