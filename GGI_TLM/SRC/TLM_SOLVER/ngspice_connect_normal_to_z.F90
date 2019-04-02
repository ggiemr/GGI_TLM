	
	if (face_update_code(face_number).NE.0) then   ! not free space

!         face_update_code points to arrays which tell us about materials, excitations and outputs
	
	  material_number=abs( face_update_code_to_material_data(face_update_code(face_number),1) )
	  if (face_update_code_to_material_data(face_update_code(face_number),1).GT.0) then
	    reverse_material=.TRUE.
	  else
	    reverse_material=.FALSE.
	  end if
	  if (material_number.NE.0) then
	    material_type=surface_material_list(material_number)%type
	  else
	    material_type=0
	  end if
	  						
	  if (material_type.NE.0) then ! not free space scatter
		  
            if (material_type.EQ.surface_material_type_SPICE) then	 
                                        
              spice_port=surface_material_list(material_number)%Spice_circuit_file_port

! connection process using the ngspice node voltage

              sign=surface_material_list(material_number)%Spice_port_sign
              
	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-x').OR.	&
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+x') ) then
  
! Transfer the TLM incident voltage pulse(s) to the ngspice circuit

                 command_string=''
                 Vspice=sign*(V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz-1))
                 if (abs(Vspice).LT.small) Vspice=0d0

! Build the command string including the spice node number                 
                 write(command_string,'(A10,I0,A3,ES16.6)')"alter vtlm",spice_port," = ",Vspice
                  
                 istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR); 
                              	      
               end if
	      
	       if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-y').OR.	&
	           (surface_material_list(material_number)%Spice_port_direction.EQ.'+y') ) then
                    
! Transfer the TLM incident voltage pulse(s) to the ngspice circuit

                 command_string=''
                 Vspice=sign*(V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz-1))
                 if (abs(Vspice).LT.small) Vspice=0d0

! Build the command string including the spice node number                 
                 write(command_string,'(A10,I0,A3,ES16.6)')"alter vtlm",spice_port," = ",Vspice
                 
                 istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR); 
                 
! ****** DEBUGGING OUTPUT TO BE REMOVED ********                
!                 write(*,*)'GGI_TLM to Spice: port',spice_port,' sign=',sign,' Vspice=',Vspice

               end if
	    
            
	     end if  ! Spice link
	    
	   end if  ! not free space scatter

	 end if  ! not free space

