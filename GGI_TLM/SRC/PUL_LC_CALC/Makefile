PUL_LC_CALC_OBJS= \
calc_D_matrix.o \
calc_Cg_matrix.o \
calc_C_matrix.o \
write_matrix.o \
calc_conductor_points.o \
calc_dielectric_points.o \
write_geometry.o \
read_geometry.o \
memory_stuff.o \
calc_phi.o \
plot_phi.o \
calc_E.o \
plot_E.o \
PUL_LC_tlm_interface.o 

default: $(PUL_LC_CALC_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
