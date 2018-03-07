default: microstrip  wire_over_ground coax two_wire

microstrip: microstrip.F90
	$(FC) $(FLAGS) -o microstrip microstrip.F90 $(LIBS)
	mv microstrip $(EXECUTABLE_DIR)/GGI_microstrip

wire_over_ground: wire_over_ground.F90
	$(FC) $(FLAGS) -o wire_over_ground wire_over_ground.F90 $(LIBS)
	mv wire_over_ground $(EXECUTABLE_DIR)/GGI_wire_over_ground

coax: coax.F90
	$(FC) $(FLAGS) -o coax coax.F90 $(LIBS)
	mv coax $(EXECUTABLE_DIR)/GGI_coax

two_wire: two_wire.F90
	$(FC) $(FLAGS) -o two_wire two_wire.F90 $(LIBS)
	mv two_wire $(EXECUTABLE_DIR)/GGI_two_wire

clean:
	rm -f *.o
	rm -f microstrip wire_over_ground coax two_wire
