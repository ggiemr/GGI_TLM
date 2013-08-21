SHELL = /bin/sh

#
# MAKEFILE NOTES:
#
# for parallel (mpi) solver only
# 'make PARALLEL_MODE=MPI'
#
# for sequential solver and all supporting codes
# 'make PARALLEL_MODE=SEQ'
#
# 'make clean' cleans everything except executables
# 'make clean_all' cleans everything including executables
#
# note that you must 'make clean' before switching from making MPI to SEQ or vice-versa

#
# May2013: these now set in the file ../configuration_data
#
#### THE FOLLOWING LINES MAY HAVE TO CHANGED ###
#
## Name of fortran compiler. Try to use one that automatically links mpi eg mpif90.
#FC=mpif90.mpich2 
#
## MPI_INCLUDE :This is the path to the file mpif.h 
#MPI_INCLUDE=/usr/include/mpich2
#
## Compilation flags, Here we compile with bounds check - slow...
##FLAGS= -D$(MODE) -fdefault-real-8 -cpp -fbounds-check -finit-real=nan -finit-integer=nan -J$(MOD_DIR) -I$(MOD_DIR) -I$(MPI_INCLUDE) 
#
## Compile without bounds check - faster...
##FLAGS= -D$(MODE) -fdefault-real-8 -cpp  -finit-real=nan -finit-integer=nan -J$(MOD_DIR) -I$(MOD_DIR) -I$(MPI_INCLUDE) 
#
## Compilation flags, Here we compile with optimisation
#FLAGS= -D$(MODE) -fdefault-real-8 -cpp -J$(MOD_DIR) -I$(MOD_DIR) -I$(MPI_INCLUDE) -O4
#
### NOTHING BELOW THIS SHOULD NEED TO BE CHANGED ###

include ../make_configuration_data

# Automatically set the paths to directories for executables, module files and libraries.

SRC_DIRECTORY:= $(shell pwd)

GGI_TLM_DIRECTORY:=$(shell dirname ${SRC_DIRECTORY})

EXECUTABLE_DIR=$(GGI_TLM_DIRECTORY)/bin

MOD_DIR=$(GGI_TLM_DIRECTORY)/GGI_TLM_MODULE_FILES

GGI_LIB_DIR=$(GGI_TLM_DIRECTORY)/GGI_TLM_LIBS

LIBS=	\
$(GGI_LIB_DIR)/lib_GGI_filter.a \
$(GGI_LIB_DIR)/lib_GGI_polynomial.a \
$(GGI_LIB_DIR)/lib_GGI_maths.a \
$(GGI_LIB_DIR)/lib_EISPACK.a 

export MPI_INCLUDE
export MOD_DIR
export FLAGS
export LIBS
export FC
export EXECUTABLE_DIR

#_____________________________________________________________

VERSION_NUMBER= 1.0.0

MODEL_BUILDER_NAME=GGI_TLM_model_builder

CABLE_MODEL_BUILDER_NAME=GGI_TLM_cable_model_builder

MODEL_CHECKS_NAME=GGI_TLM_model_checks

CABLE_MODEL_CHECKS_NAME=GGI_TLM_cable_model_checks

MATERIAL_MODEL_CHECKS_NAME=GGI_TLM_material_model_checks

FILTER_FIT_NAME=GGI_TLM_filter_fit

FILTER_FORMAT_CONVERT_NAME=GGI_TLM_filter_format_convert

SOLVER_NAME=GGI_TLM

POST_PROCESSOR_NAME=GGI_TLM_post_process

#_____________________________________________________________

CONSTANTS_OBJS=  \
CONSTANTS/constants.o

MATHS_OBJS=  \
MATHS/dmatrix.o \
MATHS/cmatrix.o \
MATHS/Legendre.o \
MATHS/bessel.o 

EISPACK_OBJS=  \
EISPACK/eispack.o 

POLYNOMIAL_OBJS=  \
POLYNOMIAL/polynomial_types.o \
POLYNOMIAL/polynomial_functions.o \
POLYNOMIAL/polynomial_operators.o \
POLYNOMIAL/root_finding.o

FILTER_OBJS=  \
FILTER/filter_types.o \
FILTER/filter_functions.o \
FILTER/filter_operators.o \
FILTER/filter_subroutines.o 

MODULE_OBJS=  \
TLM_MODULES/TLM_Modules.o \
TLM_MODULES/TLM_material_modules.o \
TLM_MODULES/TLM_excitation_modules.o \
TLM_MODULES/TLM_output_modules.o \
TLM_MODULES/TLM_mode_stir_modules.o \
TLM_MODULES/Mesh_Modules.o \
TLM_MODULES/Cable_Modules.o \
TLM_MODULES/Geometry_Modules.o \
TLM_MODULES/pul_LC_Modules.o \
TLM_MODULES/Filter_fit_modules.o \
TLM_MODULES/File_information.o \
TLM_MODULES/post_process_modules.o 

# make MODULE_OBJS available to all Makefiles 
export MODULE_OBJS

FILE_HANDLING_OBJS= \
FILE_HANDLING/file_stuff.o \
FILE_HANDLING/read_problem_name.o \
FILE_HANDLING/read_write_header_data.o 

MEMORY_HANDLING_OBJS= \
MEMORY_HANDLING/memory_stuff.o

PROCESS_INPUT_FILE_OBJS=  \
PROCESS_INPUT_FILE/reset_packet_data.o \
PROCESS_INPUT_FILE/read_input_file_geometry.o \
PROCESS_INPUT_FILE/read_input_file_solver.o \
PROCESS_INPUT_FILE/read_input_file_cables.o \
PROCESS_INPUT_FILE/read_surface_list.o \
PROCESS_INPUT_FILE/read_volume_list.o \
PROCESS_INPUT_FILE/read_line_list.o \
PROCESS_INPUT_FILE/read_point_list.o \
PROCESS_INPUT_FILE/read_mesh_dimensions_in_cells.o \
PROCESS_INPUT_FILE/read_mesh_cell_dimension.o \
PROCESS_INPUT_FILE/read_mesh_outer_boundary_dimension.o \
PROCESS_INPUT_FILE/read_simulation_time.o \
PROCESS_INPUT_FILE/read_surface_material_list.o \
PROCESS_INPUT_FILE/read_volume_material_list.o \
PROCESS_INPUT_FILE/read_excitation_function_list.o \
PROCESS_INPUT_FILE/read_excitation_point_list.o \
PROCESS_INPUT_FILE/read_excitation_surface_list.o \
PROCESS_INPUT_FILE/read_huygens_surface.o \
PROCESS_INPUT_FILE/read_far_field_surface.o \
PROCESS_INPUT_FILE/read_rcs_surface.o \
PROCESS_INPUT_FILE/read_frequency_output_surface.o \
PROCESS_INPUT_FILE/read_frequency_output_volume.o \
PROCESS_INPUT_FILE/read_frequency_domain_power_surface.o \
PROCESS_INPUT_FILE/read_output_point_list.o \
PROCESS_INPUT_FILE/read_output_surface_list.o \
PROCESS_INPUT_FILE/read_output_volume_list.o \
PROCESS_INPUT_FILE/read_output_volume_average_list.o \
PROCESS_INPUT_FILE/read_SAR_volume_list.o \
PROCESS_INPUT_FILE/read_output_time_information.o \
PROCESS_INPUT_FILE/read_outer_boundary_reflection_coefficient.o \
PROCESS_INPUT_FILE/read_field_component.o \
PROCESS_INPUT_FILE/read_centre_or_face.o \
PROCESS_INPUT_FILE/read_cable_geometry_list.o \
PROCESS_INPUT_FILE/read_cable_list.o \
PROCESS_INPUT_FILE/read_cable_output_list.o \
PROCESS_INPUT_FILE/read_cable_junction_list.o \
PROCESS_INPUT_FILE/check_cable_input_data.o \
PROCESS_INPUT_FILE/check_solver_input_data.o \
PROCESS_INPUT_FILE/read_excitation_mode_list.o \
PROCESS_INPUT_FILE/read_output_mode_list.o \
PROCESS_INPUT_FILE/read_mode_stir_surface_list.o \
PROCESS_INPUT_FILE/read_frequency_warp_flag.o \
PROCESS_INPUT_FILE/read_frequency_scale_flag.o  \
PROCESS_INPUT_FILE/read_wrapping_boundary_conditions.o

F90_GENERAL_OBJS=  \
F90_GENERAL/convert_case.o \
F90_GENERAL/strip_filename.o \
F90_GENERAL/read_write_subroutines.o \
F90_GENERAL/error_subroutines.o \
F90_GENERAL/add_integer_to_filename.o \
F90_GENERAL/write_time.o

GEOMETRY_OBJS=  \
GEOMETRY/geometry_operators.o \
GEOMETRY/transformations.o \
GEOMETRY/functions.o \
GEOMETRY/build_volume_geometry.o \
GEOMETRY/build_surface_geometry.o \
GEOMETRY/build_line_geometry.o \
GEOMETRY/build_point_geometry.o \
GEOMETRY/get_point.o \
GEOMETRY/build_surface_sphere.o \
GEOMETRY/build_surface_rectangular_block.o  \
GEOMETRY/build_surface_rectangular_block2.o \
GEOMETRY/build_surface_circle.o \
GEOMETRY/build_surface_cylinder.o \
GEOMETRY/build_surface_rectangle.o \
GEOMETRY/build_surface_xplane.o \
GEOMETRY/build_surface_yplane.o \
GEOMETRY/build_surface_zplane.o \
GEOMETRY/build_surface_quad.o \
GEOMETRY/build_surface_triangle.o \
GEOMETRY/build_surface_triangulated_surface.o \
GEOMETRY/build_surface_vtk_triangulated_surface.o \
GEOMETRY/build_volume_sphere.o \
GEOMETRY/build_volume_rectangular_block.o  \
GEOMETRY/build_volume_rectangular_block2.o \
GEOMETRY/build_volume_cylinder.o \
GEOMETRY/build_volume_tet.o \
GEOMETRY/build_volume_pyramid.o \
GEOMETRY/build_volume_pyramid_ram.o \
GEOMETRY/build_volume_tet_mesh.o \
GEOMETRY/build_line_straight_line.o \
GEOMETRY/build_line_straight_line2.o \
GEOMETRY/build_line_arc.o 

PLOTTING_OBJS=  \
PLOTTING/plot_tet_volumes.o \
PLOTTING/plot_triangulated_surfaces.o \
PLOTTING/plot_line_segments.o \
PLOTTING/plot_points.o \
PLOTTING/write_volume_list_vtk.o  \
PLOTTING/write_surface_list_vtk.o  \
PLOTTING/write_line_list_vtk.o  \
PLOTTING/write_point_list_vtk.o  \
PLOTTING/plot_mesh_volumes.o \
PLOTTING/plot_mesh_surfaces.o \
PLOTTING/plot_mesh_lines.o \
PLOTTING/plot_mesh_points.o \
PLOTTING/plot_mesh_boundary.o \
PLOTTING/write_volume_mesh_list_vtk.o \
PLOTTING/write_surface_mesh_list_vtk.o \
PLOTTING/write_line_mesh_list_vtk.o \
PLOTTING/write_point_mesh_list_vtk.o \
PLOTTING/write_tube_points_vtk.o\
PLOTTING/plot_volume_material_cells.o \
PLOTTING/plot_surface_material_faces.o

MESH_OBJS=  \
MESH/reset_mesh_parameters.o \
MESH/set_mesh_parameters.o \
MESH/general_mesh_subroutines.o \
MESH/build_volume_mesh.o \
MESH/build_surface_mesh.o \
MESH/build_line_mesh.o \
MESH/build_point_mesh.o \
MESH/mesh_triangle.o\
MESH/mesh_tet.o\
MESH/mesh_line_segment.o \
MESH/write_mesh.o \
MESH/read_mesh.o \
MESH/trim_mesh.o \
MESH/mesh_triangle_new.o \
MESH/mesh_tet_new.o \
MESH/local_grid_subroutines.o

CABLE_MODEL_OBJS=  \
CABLE_MODEL/general_cable_subroutines.o  \
CABLE_MODEL/build_cell_centre_junction_list.o  \
CABLE_MODEL/build_face_junction_list.o  \
CABLE_MODEL/create_cable_LCRG_matrices.o  \
CABLE_MODEL/create_bundle_LCRG_matrices.o  \
CABLE_MODEL/create_bundle_segments.o  \
CABLE_MODEL/get_bundle_segment_geometries.o  \
CABLE_MODEL/create_cable_meshes.o  \
CABLE_MODEL/set_bundle_excitations.o  \
CABLE_MODEL/set_bundle_outputs.o  \
CABLE_MODEL/write_cable_model_data.o \
CABLE_MODEL/read_cable_model_data.o \
CABLE_MODEL/create_bundle_cross_section_geometry.o \
CABLE_MODEL/check_and_stabilise_cable_model.o \
CABLE_MODEL/calculate_TLM_cable_matrices.o \
CABLE_MODEL/partition_cable_mesh.o \
CABLE_MODEL/cable_model_checks.o

PUL_LC_CALC_OBJS= \
PUL_LC_CALC/calc_D_matrix.o \
PUL_LC_CALC/calc_Cg_matrix.o \
PUL_LC_CALC/calc_C_matrix.o \
PUL_LC_CALC/write_matrix.o \
PUL_LC_CALC/calc_conductor_points.o \
PUL_LC_CALC/calc_dielectric_points.o \
PUL_LC_CALC/write_geometry.o \
PUL_LC_CALC/read_geometry.o \
PUL_LC_CALC/memory_stuff.o \
PUL_LC_CALC/calc_phi.o \
PUL_LC_CALC/plot_phi.o \
PUL_LC_CALC/calc_E.o \
PUL_LC_CALC/plot_E.o \
PUL_LC_CALC/PUL_LC_tlm_interface.o

FILTER_FIT_OBJS=  \
FILTER_FIT/read_filter_fit_information.o \
FILTER_FIT/read_data.o \
FILTER_FIT/filter_fit_memory_stuff.o \
FILTER_FIT/get_initial_filter_coefficients.o \
FILTER_FIT/get_testing_frequencies.o \
FILTER_FIT/optimise_filter.o \
FILTER_FIT/filter_param_conversion.o \
FILTER_FIT/stabilise_FF_input_data.o \
FILTER_FIT/stabilise_filter.o \
FILTER_FIT/stabilise_material_filter.o \
FILTER_FIT/stabilise_thin_layer_filter.o \
FILTER_FIT/stabilise_impedance_filter.o \
FILTER_FIT/test_stability.o \
FILTER_FIT/write_FF_input_data.o \
FILTER_FIT/write_filter_frequency_response.o \
FILTER_FIT/write_filter.o \
FILTER_FIT/calculate_mean_square_error.o

MATERIAL_OBJS= \
MATERIALS/TLM_volume_materials.o \
MATERIALS/TLM_surface_materials.o \
MATERIALS/output_material_frequency_response.o \
MATERIALS/Surface_material_frequency_response.o \
MATERIALS/Volume_material_frequency_response.o

EXCITATION_OBJS=  \
EXCITATION/excitation.o \
EXCITATION/excitation_points.o \
EXCITATION/excitation_surfaces.o \
EXCITATION/mode_excitation.o \
EXCITATION/huygens_surface.o 

OUTPUT_OBJS= \
OUTPUT/output.o\
OUTPUT/output_points.o \
OUTPUT/output_surfaces.o \
OUTPUT/output_volumes.o \
OUTPUT/output_volume_averages.o \
OUTPUT/output_modes.o \
OUTPUT/frequency_output_surfaces.o \
OUTPUT/frequency_output_volumes.o \
OUTPUT/frequency_domain_power_surfaces.o \
OUTPUT/far_field_surfaces.o \
OUTPUT/RCS_surfaces.o \
OUTPUT/SAR_volumes.o \
OUTPUT/cable_outputs.o

TLM_SOLVER_OBJS=  \
TLM_SOLVER/set_TLM_solver_parameters.o\
TLM_SOLVER/allocate_mesh.o \
TLM_SOLVER/set_update_codes.o \
TLM_SOLVER/run_TLM.o \
TLM_SOLVER/scatter.o \
TLM_SOLVER/connect.o \
TLM_SOLVER/TLM_parallel_pass_data.o \
TLM_SOLVER/outer_boundary.o \
TLM_SOLVER/wrap_outer_boundary.o \
TLM_SOLVER/mode_stir.o \
TLM_SOLVER/cable_bundles.o \
TLM_SOLVER/cell_centre_cable_junction.o \
TLM_SOLVER/face_cable_junction.o \
TLM_SOLVER/Cable_parallel_pass_data.o 

POST_PROCESSOR_OBJS=  \
POST_PROCESSOR/allocate_post_data.o  \
POST_PROCESSOR/find_output_files.o  \
POST_PROCESSOR/read_write_time_domain_data.o  \
POST_PROCESSOR/read_write_frequency_domain_data.o  \
POST_PROCESSOR/extract_time_domain_data.o  \
POST_PROCESSOR/plot_time_domain_data.o  \
POST_PROCESSOR/fourier_transform.o  \
POST_PROCESSOR/plot_frequency_domain_data.o \
POST_PROCESSOR/create_surface_animation.o \
POST_PROCESSOR/extract_mode.o \
POST_PROCESSOR/create_frequency_domain_animation.o \
POST_PROCESSOR/create_vector_surface_frequency_domain_animation.o \
POST_PROCESSOR/create_vector_volume_frequency_domain_animation.o \
POST_PROCESSOR/combine_frequency_domain_data.o \
POST_PROCESSOR/sum_frequency_domain_data.o \
POST_PROCESSOR/frequency_average.o \
POST_PROCESSOR/oneport.o \
POST_PROCESSOR/twoport.o \
POST_PROCESSOR/transmission_cross_section.o \
POST_PROCESSOR/IELF.o \
POST_PROCESSOR/fourier_transform_warp.o

OBJS= \
$(CONSTANTS_OBJS) \
$(MODULE_OBJS) \
$(PROCESS_INPUT_FILE_OBJS) \
$(F90_GENERAL_OBJS) \
$(GEOMETRY_OBJS) \
$(PLOTTING_OBJS) \
$(FILE_HANDLING_OBJS) \
$(MEMORY_HANDLING_OBJS) \
$(MESH_OBJS) \
$(CABLE_MODEL_OBJS) \
$(PUL_LC_CALC_OBJS) \
$(FILTER_FIT_OBJS) \
$(MATERIAL_OBJS) \
$(EXCITATION_OBJS) \
$(OUTPUT_OBJS) \
$(TLM_SOLVER_OBJS) \
$(POST_PROCESSOR_OBJS) 

# These in libraries now...
#$(MATHS_OBJS) \
#$(EISPACK_OBJS) \
#$(POLYNOMIAL_OBJS) \
#$(FILTER_OBJS) \

ifeq ($(PARALLEL_MODE),MPI)
	MODE=MPI

default: \
MAKE_constants \
MAKE_maths \
MAKE_eispack \
MAKE_polynomial \
MAKE_filter \
MAKE_modules \
MAKE_process_input_file \
MAKE_f90_general \
MAKE_geometry \
MAKE_plotting \
MAKE_file_handling \
MAKE_memory_handling \
MAKE_mesh \
MAKE_cable_model \
MAKE_pul_lc_calc \
MAKE_filter_fit \
MAKE_materials \
MAKE_excitation \
MAKE_output \
MAKE_tlm_solver \
MAKE_post_processor \
$(SOLVER_NAME) 

else
  ifeq ($(PARALLEL_MODE),SEQ)
	MODE=SEQ

default: \
MAKE_constants \
MAKE_maths \
MAKE_eispack \
MAKE_polynomial \
MAKE_filter \
MAKE_modules \
MAKE_process_input_file \
MAKE_f90_general \
MAKE_geometry \
MAKE_plotting \
MAKE_file_handling \
MAKE_memory_handling \
MAKE_mesh \
MAKE_cable_model \
MAKE_pul_lc_calc \
MAKE_filter_fit \
MAKE_materials \
MAKE_excitation \
MAKE_output \
MAKE_tlm_solver \
MAKE_post_processor \
MAKE_analytic \
$(MODEL_BUILDER_NAME) \
$(CABLE_MODEL_BUILDER_NAME) \
$(MODEL_CHECKS_NAME) \
$(CABLE_MODEL_CHECKS_NAME) \
$(MATERIAL_MODEL_CHECKS_NAME) \
$(SOLVER_NAME) \
$(FILTER_FIT_NAME) \
$(FILTER_FORMAT_CONVERT_NAME) \
$(POST_PROCESSOR_NAME) \
RUN_SCRIPT

  else
.PHONY: write_error

write_error:	
	@echo ""
	@echo "ERROR: this makefile requires an argument"
	@echo ""
	@echo " USING THE GGI_TLM MAKEFILE :"
	@echo ""
	@echo " to make the parallel (mpi) solver code only"
	@echo " make PARALLEL_MODE=MPI"
	@echo ""
	@echo " to make the sequential solver code and all other supporting codes"
	@echo " make PARALLEL_MODE=SEQ"
	@echo ""
	@echo " to clean everything except executables"
	@echo " make clean "
	@echo ""
	@echo " to clean everything including executables"
	@echo " make clean_all "
	@echo ""
	@echo " note that you must 'make clean' before switching from making MPI to SEQ or vice-versa"

  endif
endif
	
# make at top level
$(MODEL_BUILDER_NAME): $(MODEL_BUILDER_NAME).F90 $(OBJS) 
	$(FC) $(FLAGS) -o $(MODEL_BUILDER_NAME) $(MODEL_BUILDER_NAME).F90 $(OBJS) $(LIBS)
	mv $(MODEL_BUILDER_NAME) $(EXECUTABLE_DIR)
# 
$(CABLE_MODEL_BUILDER_NAME): $(CABLE_MODEL_BUILDER_NAME).F90 $(OBJS) 
	$(FC) $(FLAGS) -o $(CABLE_MODEL_BUILDER_NAME) $(CABLE_MODEL_BUILDER_NAME).F90 $(OBJS) $(LIBS) 
	mv $(CABLE_MODEL_BUILDER_NAME) $(EXECUTABLE_DIR)
# 
$(MODEL_CHECKS_NAME): $(MODEL_CHECKS_NAME).F90 $(OBJS) 
	$(FC) $(FLAGS) -o $(MODEL_CHECKS_NAME) $(MODEL_CHECKS_NAME).F90 $(OBJS) $(LIBS)
	mv $(MODEL_CHECKS_NAME) $(EXECUTABLE_DIR)	
#
$(CABLE_MODEL_CHECKS_NAME): $(CABLE_MODEL_CHECKS_NAME).F90 $(OBJS) 
	$(FC) $(FLAGS) -o $(CABLE_MODEL_CHECKS_NAME) $(CABLE_MODEL_CHECKS_NAME).F90 $(OBJS) $(LIBS) 
	mv $(CABLE_MODEL_CHECKS_NAME) $(EXECUTABLE_DIR)	
#
$(MATERIAL_MODEL_CHECKS_NAME): $(MATERIAL_MODEL_CHECKS_NAME).F90 $(OBJS) 
	$(FC) $(FLAGS) -o $(MATERIAL_MODEL_CHECKS_NAME) $(MATERIAL_MODEL_CHECKS_NAME).F90 $(OBJS) $(LIBS)
	mv $(MATERIAL_MODEL_CHECKS_NAME) $(EXECUTABLE_DIR)	
#
$(FILTER_FIT_NAME): $(FILTER_FIT_NAME).F90 $(OBJS) 
	$(FC) $(FLAGS) -o $(FILTER_FIT_NAME) $(FILTER_FIT_NAME).F90 $(OBJS) $(LIBS)
	mv $(FILTER_FIT_NAME) $(EXECUTABLE_DIR)	
#
$(FILTER_FORMAT_CONVERT_NAME): $(FILTER_FORMAT_CONVERT_NAME).F90 $(OBJS) 
	$(FC) $(FLAGS) -o $(FILTER_FORMAT_CONVERT_NAME) $(FILTER_FORMAT_CONVERT_NAME).F90 $(OBJS) $(LIBS)
	mv $(FILTER_FORMAT_CONVERT_NAME) $(EXECUTABLE_DIR)	
# 
$(SOLVER_NAME): $(SOLVER_NAME).F90 $(OBJS)
	$(FC) $(FLAGS)  -o $(SOLVER_NAME) $(SOLVER_NAME).F90 $(OBJS) $(LIBS)
	mv $(SOLVER_NAME) $(EXECUTABLE_DIR)/$(SOLVER_NAME)_$(MODE)
# 
$(POST_PROCESSOR_NAME): $(POST_PROCESSOR_NAME).F90 $(OBJS)
	$(FC) $(FLAGS) -o $(POST_PROCESSOR_NAME) $(POST_PROCESSOR_NAME).F90 $(OBJS) $(LIBS)
	mv $(POST_PROCESSOR_NAME) $(EXECUTABLE_DIR)
# 
RUN_SCRIPT: run_GGI_TLM
	cat ../run_configuration_data run_GGI_TLM > $(EXECUTABLE_DIR)/run_GGI_TLM
	chmod +x $(EXECUTABLE_DIR)/run_GGI_TLM

# make CONSTANTS
MAKE_constants:
	( cd CONSTANTS ; $(MAKE) )

# make MATHS
MAKE_maths:
	( cd MATHS ; $(MAKE) )
	ar rv $(GGI_LIB_DIR)/lib_GGI_maths.a $(MATHS_OBJS)
	
# make EISPACK
MAKE_eispack:
	( cd EISPACK ; $(MAKE) )
	ar rv $(GGI_LIB_DIR)/lib_EISPACK.a $(EISPACK_OBJS)
	
# make POLYNOMIAL
MAKE_polynomial:
	( cd POLYNOMIAL ; $(MAKE) )
	ar rv $(GGI_LIB_DIR)/lib_GGI_polynomial.a $(POLYNOMIAL_OBJS)

# make FILTER
MAKE_filter:
	( cd FILTER ; $(MAKE) )
	ar rv $(GGI_LIB_DIR)/lib_GGI_filter.a $(FILTER_OBJS)

# make MODULES
MAKE_modules:
	( cd TLM_MODULES ; $(MAKE) )

# make PROCESS_INPUT_FILE
MAKE_process_input_file:
	( cd PROCESS_INPUT_FILE ; $(MAKE) )

# make F90_GENERAL
MAKE_f90_general:
	( cd F90_GENERAL ; $(MAKE) )

# make GEOMETRY
MAKE_geometry:
	( cd GEOMETRY ; $(MAKE) )

# make PLOTTING
MAKE_plotting:
	( cd PLOTTING ; $(MAKE) )

# make FILE_HANDLING
MAKE_file_handling:
	( cd FILE_HANDLING ; $(MAKE) )

# make MEMORY_HANDLING
MAKE_memory_handling:
	( cd MEMORY_HANDLING ; $(MAKE) )

# make MESH
MAKE_mesh:
	( cd MESH ; $(MAKE) )
	
# make CABLE_MODEL
MAKE_cable_model:
	( cd CABLE_MODEL ; $(MAKE) )
	
# make PUL_LC
MAKE_pul_lc_calc:
	( cd PUL_LC_CALC ; $(MAKE) )

# make FILTER_FIT
MAKE_filter_fit:
	( cd FILTER_FIT ; $(MAKE) )
	
# make MATERIALS
MAKE_materials:
	( cd MATERIALS ; $(MAKE) )

# make EXCITATION
MAKE_excitation:
	( cd EXCITATION ; $(MAKE) )

# make OUTPUT
MAKE_output:
	( cd OUTPUT ; $(MAKE) )

# make TLM_SOLVER
MAKE_tlm_solver:
	( cd TLM_SOLVER ; $(MAKE) )
	
# make POST_PROCESSOR
MAKE_post_processor:
	( cd POST_PROCESSOR ; $(MAKE) )
	
# make ANALYTIC
MAKE_analytic:
	( cd ANALYTIC ; $(MAKE) )

clean:
	( rm -f *.o )
	( rm -f $(MOD_DIR)/*.mod )
	( rm -f $(GGI_LIB_DIR)/*.a )
	( cd CONSTANTS ; make clean )
	( cd MATHS ; make clean )
	( cd EISPACK ; make clean )
	( cd POLYNOMIAL ; make clean )
	( cd FILTER ; make clean )
	( cd TLM_MODULES ; make clean )
	( cd PROCESS_INPUT_FILE ; make clean )
	( cd F90_GENERAL ; make clean )
	( cd GEOMETRY ; make clean )
	( cd PLOTTING ; make clean )
	( cd FILE_HANDLING ; make clean )
	( cd MEMORY_HANDLING ; make clean )
	( cd MESH ; make clean )
	( cd CABLE_MODEL ; make clean )
	( cd MATERIALS ; make clean )
	( cd PUL_LC_CALC ; make clean )
	( cd FILTER_FIT ; make clean )
	( cd EXCITATION ; make clean )
	( cd OUTPUT ; make clean )
	( cd TLM_SOLVER ; make clean )
	( cd POST_PROCESSOR ; make clean )
	( cd ANALYTIC ; make clean )

clean_all:
	make clean
	( rm -f $(EXECUTABLE_DIR)/*GGI_* )



