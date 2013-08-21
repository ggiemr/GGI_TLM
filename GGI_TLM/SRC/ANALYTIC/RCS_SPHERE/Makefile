default: RCS_sphere near_field_convolution

OBJS=	RCS_sphere.o tests.o rcs.o rcs_ibc.o rcs_coated.o rcs_coated_ibc.o rcs_dielectric_ibc.o
#  
RCS_sphere: $(OBJS)
	$(FC) $(FLAGS) -o RCS_sphere $(OBJS) $(LIBS) 
	mv RCS_sphere $(EXECUTABLE_DIR)/GGI_RCS_sphere
#  
near_field_convolution: near_field_convolution.o
	$(FC) $(FLAGS) -o near_field_convolution near_field_convolution.o $(LIBS) 
	mv near_field_convolution $(EXECUTABLE_DIR)/GGI_near_field_convolution
#

%.o: %.F90
	$(FC) $(FLAGS) -c $<
#
clean:
	rm -f *.o
	
