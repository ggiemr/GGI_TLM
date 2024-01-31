FC=f95
default: wire_over_lossy_ground

wire_over_lossy_ground: wire_over_lossy_ground.F90 calc_error_Wait.F90 mod_zbes.o
	$(FC) $(FLAGS) -o  wire_over_lossy_ground wire_over_lossy_ground.F90 mod_zbes.o $(LIBS)
	mv wire_over_lossy_ground $(EXECUTABLE_DIR)/GGI_wire_over_lossy_ground
	
mod_zbes.o: mod_zbes.f90
	$(FC) $(FLAGS) -c mod_zbes.f90

clean:
	rm -f *.o *.mod *.out
	rm -f wire_over_lossy_ground
