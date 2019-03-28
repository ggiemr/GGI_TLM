	
	if (face_update_code(face_number).eq.0) then   ! free space update, no excitation, no output
	
	  Vx=V(Vx_zmin,cx,cy,cz)+ V(Vx_zmax,cx,cy,cz-1)
	  V(Vx_zmin,cx,cy,cz)  =  Vx-V(Vx_zmin,cx,cy,cz)
	  V(Vx_zmax,cx,cy,cz-1)=  Vx-V(Vx_zmax,cx,cy,cz-1)
       
	  Vy=V(Vy_zmin,cx,cy,cz)+ V(Vy_zmax,cx,cy,cz-1)
	  V(Vy_zmin,cx,cy,cz)  =  Vy-V(Vy_zmin,cx,cy,cz)
	  V(Vy_zmax,cx,cy,cz-1)=  Vy-V(Vy_zmax,cx,cy,cz-1)	
	  
	else
!         face_update_code points to arrays which tell us about material number and output_number
	
	  material_number=abs(face_update_code_to_material_data(face_update_code(face_number),1))
	  	  
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
	  
	  cable_face_junction_number=face_update_code_to_cable_number(face_update_code(face_number))
	  output_number=face_update_code_to_output_number(face_update_code(face_number))
	  excitation_number=abs( face_update_code_to_excitation_number(face_update_code(face_number)) )
			
	  if (output_number.NE.0) then	  
	  
! get the incident voltage pulses required for output field calculation	  
	    Vx_negative_z1=V(Vx_zmin,cx,cy,cz)
	    Vy_negative_z1=V(Vy_zmin,cx,cy,cz)
	    Vx_positive_z2=V(Vx_zmax,cx,cy,cz-1)
	    Vy_positive_z2=V(Vy_zmax,cx,cy,cz-1)
	    
	  end if		
	
	  if (material_type.EQ.0) then ! no material, simple connect
	
	    Vx_min=V(Vx_zmin,cx,cy,cz)+ V(Vx_zmax,cx,cy,cz-1)
	    Vx_max=Vx_min
       
	    Vy_max=V(Vy_zmin,cx,cy,cz)+ V(Vy_zmax,cx,cy,cz-1)
	    Vy_min=Vy_max
	  
	  else ! material_type.NE.0
	  
	    if (material_type.EQ.surface_material_type_PEC) then
	  
	      Vx_min=0D0
	      Vx_max=0D0
         
	      Vy_max=0D0
	      Vy_min=0D0
	  
	    else if (material_type.EQ.surface_material_type_PMC) then
	  
	      Vx_min=2D0*V(Vx_zmin,cx,cy,cz)
	      Vx_max=2D0*V(Vx_zmax,cx,cy,cz-1)
       
	      Vy_min=2D0*V(Vy_zmin,cx,cy,cz)
	      Vy_max=2D0*V(Vy_zmax,cx,cy,cz-1)
	  
	    else if (material_type.EQ.surface_material_type_FREE_SPACE) then	 
	  
	      Vx_min=V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz-1)
	      Vx_max=Vx_min
       
	      Vy_min=V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz-1)
	      Vy_max=Vy_min
	  
	    else if (material_type.EQ.surface_material_type_DIODE) then	 
	    
	      Is =surface_material_list(material_number)%Diode_Is
	      nVt=surface_material_list(material_number)%Diode_nVt
	      Rs=surface_material_list(material_number)%Diode_Rs
              sign=surface_material_list(material_number)%Diode_sign  
	      
	      if ((surface_material_list(material_number)%Diode_direction.EQ.'-x').OR.	&
	          (surface_material_list(material_number)%Diode_direction.EQ.'+x') ) then

! diode junction capacitance TLM model update	    
                diode_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)
                CALL timeshift_Zfilter(Diode_Cj_filter_data(diode_filter_number))
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       0d0)

                Zc=surface_material_list(material_number)%Diode_Cj_f
		Vc=Diode_Cj_filter_data(diode_filter_number)%f
	      
	        CALL diode_calc(V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz-1),Z0/2d0,Is,nVt,Rs,Zc,Vc,Vd,Id,Ic,sign)
		Vx_min=Vd
		Vx_max=Vd
		
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       Ic)
	      
	      else ! diode is not is this direction so do free space update for x polarisation
	      
	        Vx_min=V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz-1)
	        Vx_max=Vx_min
       
              end if
	      
	      if ((surface_material_list(material_number)%Diode_direction.EQ.'-y').OR.	&
	          (surface_material_list(material_number)%Diode_direction.EQ.'+y') ) then

! diode junction capacitance TLM model update	    
                diode_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)
                CALL timeshift_Zfilter(Diode_Cj_filter_data(diode_filter_number))
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       0d0)

                Zc=surface_material_list(material_number)%Diode_Cj_f
		Vc=Diode_Cj_filter_data(diode_filter_number)%f
	      
	        CALL diode_calc(V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz-1),Z0/2d0,Is,nVt,Rs,Vd,Id,1)
		Vy_min=Vd
		Vy_max=Vd
		
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       Ic)
	      	      
	      else ! diode is not is this direction so do free space update for y polarisation
	      
	        Vy_min=V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz-1)
	        Vy_max=Vy_min
       
              end if	      
	  
	    else if (material_type.EQ.surface_material_type_SPICE) then	 
                  
              spice_node1=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
              spice_node2=surface_material_list(material_number)%Spice_circuit_file_nodes(2)

! connection process using the ngspice node voltage
              sign=surface_material_list(material_number)%Spice_port_sign

	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-x').OR.	&
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+x') ) then
             
! get the ngspice node voltage. The voltage to use is found in the list ngspice_node_to_V_ngspice_array_list(100)
                opnode1=ngspice_node_to_V_ngspice_array_list(spice_node1)
                opnode2=ngspice_node_to_V_ngspice_array_list(spice_node2)
                                
                if ( (opnode1.LT.0).OR.(opnode1.GT.100) ) then
                  write(*,*)'Ngspice output node1 is out of range (0-100)',opnode1
                  STOP
                end if
                if ( (opnode2.LT.0).OR.(opnode2.GT.100) ) then
                  write(*,*)'Ngspice output node2 is out of range (0-100)',opnode2
                  STOP
                end if
                
                if (opnode1.NE.0) then
                  Vspice1=sign*V_ngspice_array_F90(opnode1) 
                else
                  Vspice1=0d0
                end if 
                
                if (opnode2.NE.0) then
                  Vspice2=V_ngspice_array_F90(opnode2) 
                else
                  Vspice2=0d0
                end if 
                Vspice=sign*(Vspice1-Vspice2)    
                
		Vx_min=Vspice
		Vx_max=Vspice
	      
	      else ! spice circuit port is not is this direction so do free space update for x polarisation
	      
	        Vx_min=V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz-1)
	        Vx_max=Vx_min

              end if
	      
	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-y').OR.	&
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+y') ) then
                              
! get the ngspice node voltage. The voltage to use is found in the list ngspice_node_to_V_ngspice_array_list(100)
                opnode1=ngspice_node_to_V_ngspice_array_list(spice_node1)
                opnode2=ngspice_node_to_V_ngspice_array_list(spice_node2)
                
                if ( (opnode1.LT.0).OR.(opnode1.GT.100) ) then
                  write(*,*)'Ngspice output node1 is out of range (0-100)',opnode1
                  STOP
                end if
                if ( (opnode2.LT.0).OR.(opnode2.GT.100) ) then
                  write(*,*)'Ngspice output node2 is out of range (0-100)',opnode2
                  STOP
                end if
                
                if (opnode1.NE.0) then
                  Vspice1=sign*V_ngspice_array_F90(opnode1) 
                else
                  Vspice1=0d0
                end if 
                
                if (opnode2.NE.0) then
                  Vspice2=V_ngspice_array_F90(opnode2) 
                else
                  Vspice2=0d0
                end if 
                Vspice=sign*(Vspice1-Vspice2)    
                                 
                Vy_min=Vspice
		Vy_max=Vspice

! ****** DEBUGGING OUTPUT TO BE REMOVED ********                
!                write(*,*)'Spice to GGI_TLM: port',surface_material_list(material_number)%Spice_circuit_file_port, &
!                          ' nodes',opnode1,opnode2,' Vspice=',Vspice

	      else ! spice circuit port is not is this direction so do free space update for y polarisation
	      
	        Vy_min=V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz-1)
	        Vy_max=Vy_min
                
              end if
	     
	    else if ( (material_type.EQ.surface_material_type_DISPERSIVE).OR.			&
	              (material_type.EQ.surface_material_type_ANISOTROPIC_DISPERSIVE) ) then
	    
	      surface_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)

! Vx polarisation 	      
              pol=1
              call surface_material_update(V(Vx_zmax,cx,cy,cz-1),Z0,V(Vx_zmin,cx,cy,cz),Z0,	&
	                                   Vx_max,Vx_min,material_number,pol,reverse_material,surface_filter_number)	
	      	      
! Vy polarisation 	      
              pol=2
              call surface_material_update(V(Vy_zmax,cx,cy,cz-1),Z0,V(Vy_zmin,cx,cy,cz),Z0,	&
	                                   Vy_max,Vy_min,material_number,pol,reverse_material,surface_filter_number+1)	
	  
	    else if (material_type.EQ.surface_material_type_SPICE) then	 
	  
	    end if
	    
	  end if

! Excitation		
          if (excitation_number.ne.0) then
	  
	    excitation_face_number=excitation_face_number+1
	    field_min(1:6)=face_excitation_field(excitation_face_number,1,1:6)
	    field_max(1:6)=face_excitation_field(excitation_face_number,2,1:6)
	    
	    hard_source_factor_min(1:6)=face_excitation_type(excitation_face_number,1,1:6)
	    hard_source_factor_max(1:6)=face_excitation_type(excitation_face_number,1,1:6)

	    Js_min(1:3)=face_excitation_field(excitation_face_number,1,7:9)
	    Js_max(1:3)=face_excitation_field(excitation_face_number,2,7:9)
	    Ms_min(1:3)=face_excitation_field(excitation_face_number,1,10:12)
	    Ms_max(1:3)=face_excitation_field(excitation_face_number,2,10:12)
	    
            Vx_min=Vx_min*hard_source_factor_min(Ex)-field_min(Ex)*dl-field_min(Hy)*Z0*dl+Js_min(1)*dl*Z0/2d0+Ms_min(2)*dl/2d0
	    Vx_max=Vx_max*hard_source_factor_max(Ex)-field_max(Ex)*dl+field_max(Hy)*Z0*dl+Js_max(1)*dl*Z0/2d0-Ms_max(2)*dl/2d0
	    
            Vy_min=Vy_min*hard_source_factor_min(Ey)-field_min(Ey)*dl+field_min(Hx)*Z0*dl+Js_min(2)*dl*Z0/2d0-Ms_min(1)*dl/2d0
	    Vy_max=Vy_max*hard_source_factor_max(Ey)-field_max(Ey)*dl-field_max(Hx)*Z0*dl+Js_max(2)*dl*Z0/2d0+Ms_max(1)*dl/2d0
	  
	  end if
					       
! complete the connection					       
          V(Vx_zmax,cx,cy,cz-1)  =  Vx_max-V(Vx_zmax,cx,cy,cz-1)
	  V(Vx_zmin,cx,cy,cz  )  =  Vx_min-V(Vx_zmin,cx,cy,cz  )
	    				   
          V(Vy_zmax,cx,cy,cz-1)  =  Vy_max-V(Vy_zmax,cx,cy,cz-1)
	  V(Vy_zmin,cx,cy,cz  )  =  Vy_min-V(Vy_zmin,cx,cy,cz  )

! Cable connect		
          if (cable_face_junction_number.ne.0) then
	    
	    call Face_Cable_Junction( cable_face_junction_number )
	  
	  end if

! Output		
          if (output_number.ne.0) then
	  
! get the scattered voltage pulses required for output field calculation	  
	    Vx_positive_z1=V(Vx_zmin,cx,cy,cz)
	    Vy_positive_z1=V(Vy_zmin,cx,cy,cz)
	    Vx_negative_z2=V(Vx_zmax,cx,cy,cz-1)
	    Vy_negative_z2=V(Vy_zmax,cx,cy,cz-1)
	    
	    output_face_number=output_face_number+1
	    
! field on side 1 (zmin)
! E=-V/dl	    
	    field(Ex)=-(Vx_positive_z1+Vx_negative_z1)/dl
	    field(Ey)=-(Vy_positive_z1+Vy_negative_z1)/dl
	    field(Ez)=0d0
! H= I/dl	    
	    field(Hx)= (Vy_positive_z1-Vy_negative_z1)/(Z0*dl)
	    field(Hy)=-(Vx_positive_z1-Vx_negative_z1)/(Z0*dl)
	    field(Hz)=0d0 
	  
! Save output on this face
	    
            face_output_field(output_face_number,1,1:6)=field(1:6)	  	    
	    
! field on side 2 (zmax)
! E=-V/dl	    
	    field(Ex)=-(Vx_positive_z2+Vx_negative_z2)/dl
	    field(Ey)=-(Vy_positive_z2+Vy_negative_z2)/dl
	    field(Ez)=0d0
! H= I/dl	    
	    field(Hx)= (Vy_positive_z2-Vy_negative_z2)/(Z0*dl)
	    field(Hy)=-(Vx_positive_z2-Vx_negative_z2)/(Z0*dl)
	    field(Hz)=0d0 
	  
! Save output on this face
	    
            face_output_field(output_face_number,2,1:6)=field(1:6)	
	  
	  end if
	
	end if
