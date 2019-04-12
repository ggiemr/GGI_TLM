
TLM_SOLVER_OBJS=  \
set_TLM_solver_parameters.o \
allocate_mesh.o \
set_update_codes.o \
run_TLM.o \
run_TLM_periodic_BC.o \
scatter.o \
ngspice_init_connect.o \
connect.o \
Cable_parallel_pass_data.o \
outer_boundary.o \
wrap_outer_boundary.o \
mode_stir.o \
cable_bundles.o \
cell_centre_cable_junction.o \
face_cable_junction.o \
TLM_parallel_pass_data.o 

default: $(TLM_SOLVER_OBJS) 

%.o: %.F90 connect_normal_to_x.F90  connect_normal_to_y.F90  connect_normal_to_z.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o

