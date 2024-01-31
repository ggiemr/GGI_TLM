default: \
rectangular_waveguide \
cavity \
mlayer \
skin_depth \
rcs_cylinder \
rcs_sphere \
cable_loss \
transmission_line \
wire_over_lossy_ground

clean: \
clean_rectangular_waveguide \
clean_cavity \
clean_mlayer \
clean_skin_depth \
clean_rcs_cylinder \
clean_rcs_sphere \
clean_cable_loss \
clean_transmission_line \
clean_wire_over_lossy_ground

rectangular_waveguide:
	( cd RECTANGULAR_WAVEGUIDE ; $(MAKE) )
cavity:
	( cd CAVITY ; $(MAKE) )
mlayer:
	( cd MLAYER ; $(MAKE) )
skin_depth:
	( cd SKIN_DEPTH ; $(MAKE) )
rcs_cylinder:
	( cd RCS_CYLINDER ; $(MAKE) )
rcs_sphere:
	( cd RCS_SPHERE ; $(MAKE) )
cable_loss:
	( cd CABLE_LOSS ; $(MAKE) )
transmission_line:
	( cd TRANSMISSION_LINE ; $(MAKE) )
wire_over_lossy_ground:
	( cd WIRE_OVER_LOSSY_GROUND ; $(MAKE) )
        
#
clean_rectangular_waveguide:
	( cd RECTANGULAR_WAVEGUIDE ; make clean )
clean_cavity:
	( cd CAVITY ; make clean )
clean_mlayer:
	( cd MLAYER ; make clean )
clean_skin_depth:
	( cd SKIN_DEPTH ; make clean )
clean_rcs_cylinder:
	( cd RCS_CYLINDER ; make clean )
clean_rcs_sphere:
	( cd RCS_SPHERE ; make clean )
clean_cable_loss:
	( cd CABLE_LOSS ; make clean )
clean_transmission_line:
	( cd TRANSMISSION_LINE ; make clean )
clean_wire_over_lossy_ground:
	( cd WIRE_OVER_LOSSY_GROUND ; make clean )
