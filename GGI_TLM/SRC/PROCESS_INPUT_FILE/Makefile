
PROCESS_INPUT_FILE_OBJS=  \
reset_packet_data.o \
read_input_file_geometry.o \
read_input_file_solver.o \
read_input_file_cables.o \
read_surface_list.o \
read_volume_list.o \
read_line_list.o \
read_point_list.o \
read_mesh_dimensions_in_cells.o \
read_mesh_cell_dimension.o \
read_mesh_outer_boundary_dimension.o \
read_simulation_time.o \
read_surface_material_list.o \
read_volume_material_list.o \
read_excitation_function_list.o \
read_excitation_point_list.o \
read_excitation_surface_list.o \
read_huygens_surface.o \
read_far_field_surface.o \
read_rcs_surface.o \
read_frequency_output_surface.o \
read_frequency_output_volume.o \
read_frequency_domain_power_surface.o \
read_output_point_list.o \
read_output_surface_list.o \
read_output_volume_list.o \
read_output_volume_average_list.o \
read_SAR_volume_list.o \
read_output_time_information.o \
read_outer_boundary_reflection_coefficient.o \
read_field_component.o \
read_centre_or_face.o \
read_cable_geometry_list.o \
read_cable_list.o \
read_cable_output_list.o \
read_cable_junction_list.o \
check_cable_input_data.o \
check_solver_input_data.o  \
read_excitation_mode_list.o \
read_output_mode_list.o \
read_mode_stir_surface_list.o \
read_frequency_warp_flag.o \
read_frequency_scale_flag.o \
read_wrapping_boundary_conditions.o

default: $(PROCESS_INPUT_FILE_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
