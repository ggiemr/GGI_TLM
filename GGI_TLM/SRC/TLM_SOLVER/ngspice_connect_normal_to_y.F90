	
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
                 Vspice=sign*(V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy-1,cz))
                 if (abs(Vspice).LT.small) Vspice=0d0

! Apply a low pass filter to the voltage data going from TLM to Ngspice 
                 V_tlm_to_ngspice_in(spice_port,1)=Vspice      
! Apply first order LPF  fout(t)=(fin(t-1)+fin(t)-k2*fout(t-1))/k1
                 V_tlm_to_ngspice_out(spice_port,1)=(V_tlm_to_ngspice_in(spice_port,2)+V_tlm_to_ngspice_in(spice_port,1)     &
                                                                                  -LPF_k2*V_tlm_to_ngspice_out(spice_port,2))/LPF_k1     
                 Vspice=V_tlm_to_ngspice_out(spice_port,1)      
! timeshift voltage pulses
                 V_tlm_to_ngspice_in(spice_port,2)=V_tlm_to_ngspice_in(spice_port,1)
                 V_tlm_to_ngspice_out(spice_port,2)=V_tlm_to_ngspice_out(spice_port,1)

! Build the command string including the spice node number                 
                 write(command_string,'(A10,I0,A3,ES16.6)')"alter vtlm",spice_port," = ",Vspice
                 
                 istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR); 
                              	      
               end if
                                   
	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-z').OR.  &
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+z') ) then
  
! Transfer the TLM incident voltage pulse(s) to the ngspice circuit

                 command_string=''
                 Vspice=sign*(V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy-1,cz))
                 if (abs(Vspice).LT.small) Vspice=0d0

! Apply a low pass filter to the voltage data going from TLM to Ngspice 
                 V_tlm_to_ngspice_in(spice_port,1)=Vspice      
! Apply first order LPF  fout(t)=(fin(t-1)+fin(t)-k2*fout(t-1))/k1
                 V_tlm_to_ngspice_out(spice_port,1)=(V_tlm_to_ngspice_in(spice_port,2)+V_tlm_to_ngspice_in(spice_port,1)     &
                                                                                  -LPF_k2*V_tlm_to_ngspice_out(spice_port,2))/LPF_k1     
                 Vspice=V_tlm_to_ngspice_out(spice_port,1)      
! timeshift voltage pulses
                 V_tlm_to_ngspice_in(spice_port,2)=V_tlm_to_ngspice_in(spice_port,1)
                 V_tlm_to_ngspice_out(spice_port,2)=V_tlm_to_ngspice_out(spice_port,1)

! Build the command string including the spice node number                 
                 write(command_string,'(A10,I0,A3,ES16.6)')"alter vtlm",spice_port," = ",Vspice
                 
                 istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR); 
                              	      
               end if
	                
	    end if  ! Spice link
	    
	  end if  ! not free space scatter

	end if  ! not free space
