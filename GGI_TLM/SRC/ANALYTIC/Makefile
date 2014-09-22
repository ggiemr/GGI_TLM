default: \
cavity \
mlayer \
skin_depth \
rcs_cylinder \
rcs_sphere \
cable_loss \
transmission_line

clean: \
clean_cavity \
clean_mlayer \
clean_skin_depth \
clean_rcs_cylinder \
clean_rcs_sphere \
clean_cable_loss \
clean_transmission_line

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
#
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
