	
	if (face_update_code(face_number).eq.0) then   ! free space update, no excitation, no output
	
          Vy=V(Vy_xmin,cx,cy,cz)  + V(Vy_xmax,cx-1,cy,cz)
	  V(Vy_xmin,cx,cy,cz)    =  Vy-V(Vy_xmin,cx,cy,cz)
          V(Vy_xmax,cx-1,cy,cz)  =  Vy-V(Vy_xmax,cx-1,cy,cz)
	
          Vz=V(Vz_xmin,cx,cy,cz)  +  V(Vz_xmax,cx-1,cy,cz)
	  V(Vz_xmin,cx,cy,cz   )  =  Vz-V(Vz_xmin,cx,cy,cz  )
          V(Vz_xmax,cx-1,cy,cz)   =  Vz-V(Vz_xmax,cx-1,cy,cz)
	  
	else
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
	  
	  cable_face_junction_number=face_update_code_to_cable_number(face_update_code(face_number))
	  output_number=face_update_code_to_output_number(face_update_code(face_number))
	  excitation_number=abs( face_update_code_to_excitation_number(face_update_code(face_number)) )
	  	  			
	  if (output_number.NE.0) then	  
	  
! get the incident voltage pulses required for output field calculation	  
	    Vy_negative_x1=V(Vy_xmin,cx,cy,cz)
	    Vz_negative_x1=V(Vz_xmin,cx,cy,cz)
	    Vy_positive_x2=V(Vy_xmax,cx-1,cy,cz)
	    Vz_positive_x2=V(Vz_xmax,cx-1,cy,cz)
	    
	  end if		
						
	  if (material_type.EQ.0) then ! no material, simple connect
	
            Vy_min=V(Vy_xmin,cx,cy,cz)  + V(Vy_xmax,cx-1,cy,cz)
	    Vy_max=Vy_min
	    
            Vz_min=V(Vz_xmin,cx,cy,cz)  +  V(Vz_xmax,cx-1,cy,cz)
            Vz_max=Vz_min
	  
	  else ! material_type.NE.0
	  
	    if (material_type.EQ.surface_material_type_PEC) then
	
              Vy_min=0d0
	      Vy_max=0d0
	    
              Vz_min=0d0
              Vz_max=0d0
	  
	    else if (material_type.EQ.surface_material_type_PMC) then	  
	
              Vy_min=2D0*V(Vy_xmin,cx,cy,cz)  
	      Vy_max=2D0*V(Vy_xmax,cx-1,cy,cz)
	    
              Vz_min=2D0*V(Vz_xmin,cx,cy,cz) 
	      Vz_max=2D0*V(Vz_xmax,cx-1,cy,cz)
	  
	    else if (material_type.EQ.surface_material_type_FREE_SPACE) then	  
	
              Vy_min=V(Vy_xmin,cx,cy,cz)  + V(Vy_xmax,cx-1,cy,cz)
	      Vy_max=Vy_min
	    
              Vz_min=V(Vz_xmin,cx,cy,cz)  +  V(Vz_xmax,cx-1,cy,cz)
	      Vz_max=Vz_min
	  
	    else if (material_type.EQ.surface_material_type_DIODE) then	 
	    
	      Is =surface_material_list(material_number)%Diode_Is
	      nVt=surface_material_list(material_number)%Diode_nVt
	      Rs=surface_material_list(material_number)%Diode_Rs
              sign=surface_material_list(material_number)%Diode_sign  
	      
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
  
	        CALL diode_calc(V(Vy_xmin,cx,cy,cz) + V(Vy_xmax,cx-1,cy,cz),Z0/2d0,Is,nVt,Rs,Zc,Vc,Vd,Id,Ic,sign)
		Vy_min=Vd
		Vy_max=Vd
		
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       Ic)
	      
	      else ! diode is not is this direction so do free space update for y polarisation
	      
                Vy_min=V(Vy_xmin,cx,cy,cz)  + V(Vy_xmax,cx-1,cy,cz)
	        Vy_max=Vy_min
       
              end if	      
	      
	      if ((surface_material_list(material_number)%Diode_direction.EQ.'-z').OR.	&
	          (surface_material_list(material_number)%Diode_direction.EQ.'+z') ) then

! diode junction capacitance TLM model update	    
                diode_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)
                CALL timeshift_Zfilter(Diode_Cj_filter_data(diode_filter_number))
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       0d0)

                Zc=surface_material_list(material_number)%Diode_Cj_f
		Vc=Diode_Cj_filter_data(diode_filter_number)%f
	      
	        CALL diode_calc(V(Vz_xmin,cx,cy,cz)  +  V(Vz_xmax,cx-1,cy,cz),Z0/2d0,Is,nVt,Rs,Zc,Vc,Vd,Id,Ic,sign)
		Vz_min=Vd
		Vz_max=Vd
		
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       Ic)
	      
	      else ! diode is not is this direction so do free space update for z polarisation
	      
                Vz_min=V(Vz_xmin,cx,cy,cz)  +  V(Vz_xmax,cx-1,cy,cz)
	        Vz_max=Vz_min
       
              end if
	  
	    else if (material_type.EQ.surface_material_type_SPICE) then	 

#if defined(INCLUDE_NGSPICE)    

              spice_node1=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
              spice_node2=surface_material_list(material_number)%Spice_circuit_file_nodes(2)

! connection process using the ngspice node voltage
              sign=surface_material_list(material_number)%Spice_port_sign

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
                  Vspice1=V_ngspice_to_tlm_out(opnode1,1)
                else
                  Vspice1=0d0
                end if 
                
                if (opnode2.NE.0) then
                  Vspice2=V_ngspice_to_tlm_out(opnode2,1)
                else
                  Vspice2=0d0
                end if 
                Vspice=sign*(Vspice1-Vspice2)    
                
		Vy_min=Vspice
		Vy_max=Vspice
	      
	      else ! spice circuit port is not is this direction so do free space update for x polarisation
	      
                Vy_min=V(Vy_xmin,cx,cy,cz)  + V(Vy_xmax,cx-1,cy,cz)
	        Vy_max=Vy_min

              end if
	      
	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-z').OR.	&
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+z') ) then
                              
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
                  Vspice1=V_ngspice_to_tlm_out(opnode1,1)
                else
                  Vspice1=0d0
                end if 
                
                if (opnode2.NE.0) then
                  Vspice2=V_ngspice_to_tlm_out(opnode2,1)
                else
                  Vspice2=0d0
                end if 
                Vspice=sign*(Vspice1-Vspice2)    
                
                Vz_min=Vspice
		Vz_max=Vspice

	      else ! spice circuit port is not is this direction so do free space update for y polarisation
	      
                Vz_min=V(Vz_xmin,cx,cy,cz)  +  V(Vz_xmax,cx-1,cy,cz)
	        Vz_max=Vz_min
                
              end if

#endif
	  
	    else if ( (material_type.EQ.surface_material_type_DISPERSIVE).OR.			&
	              (material_type.EQ.surface_material_type_ANISOTROPIC_DISPERSIVE) ) then
	    
! Vy polarisation 	      
	      surface_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)
	      pol=2
              call surface_material_update(V(Vy_xmax,cx-1,cy,cz),Z0,V(Vy_xmin,cx,cy,cz),Z0,	&
	                                   Vy_max,Vy_min,material_number,pol,reverse_material,surface_filter_number)	
	      
! Vz polarisation 	      
	      pol=3
              call surface_material_update(V(Vz_xmax,cx-1,cy,cz),Z0,V(Vz_xmin,cx,cy,cz),Z0,	&
	                                   Vz_max,Vz_min,material_number,pol,reverse_material,surface_filter_number+1)	            
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
	    
            Vy_min=Vy_min*hard_source_factor_min(Ey)-field_min(Ey)*dl-field_min(Hz)*Z0*dl+Js_min(2)*dl*Z0/2d0+Ms_min(3)*dl/2d0
	    Vy_max=Vy_max*hard_source_factor_max(Ey)-field_max(Ey)*dl+field_max(Hz)*Z0*dl+Js_max(2)*dl*Z0/2d0-Ms_max(3)*dl/2d0
	    
            Vz_min=Vz_min*hard_source_factor_min(Ez)-field_min(Ez)*dl+field_min(Hy)*Z0*dl+Js_min(3)*dl*Z0/2d0-Ms_min(2)*dl/2d0
	    Vz_max=Vz_max*hard_source_factor_max(Ez)-field_max(Ez)*dl-field_max(Hy)*Z0*dl+Js_max(3)*dl*Z0/2d0+Ms_max(2)*dl/2d0
	    
	  end if

! complete the connection					       
          V(Vy_xmax,cx-1,cy,cz)  =  Vy_max-V(Vy_xmax,cx-1,cy,cz)
	  V(Vy_xmin,cx,cy,cz)    =  Vy_min-V(Vy_xmin,cx,cy,cz)
					       
          V(Vz_xmax,cx-1,cy,cz)  =  Vz_max-V(Vz_xmax,cx-1,cy,cz)
	  V(Vz_xmin,cx,cy,cz)    =  Vz_min-V(Vz_xmin,cx,cy,cz)	      

! Cable connect		
          if (cable_face_junction_number.ne.0) then
	    
	    CALL face_cable_junction( cable_face_junction_number )
	  
	  end if

! Output					
	  if (output_number.NE.0) then
	  
! get the scattered voltage pulses required for output field calculation	  
	    Vy_positive_x1=V(Vy_xmin,cx,cy,cz)
	    Vz_positive_x1=V(Vz_xmin,cx,cy,cz)
	    Vy_negative_x2=V(Vy_xmax,cx-1,cy,cz)
	    Vz_negative_x2=V(Vz_xmax,cx-1,cy,cz)
	    
	    output_face_number=output_face_number+1

! field on side 1, (xmin)
! E=-V/dl	    
	    field(Ex)=0d0
	    field(Ey)=-(Vy_positive_x1+Vy_negative_x1)/dl
	    field(Ez)=-(Vz_positive_x1+Vz_negative_x1)/dl
! H= I/dl	    
	    field(Hx)=0d0
	    field(Hy)= (Vz_positive_x1-Vz_negative_x1)/(Z0*dl)
	    field(Hz)=-(Vy_positive_x1-Vy_negative_x1)/(Z0*dl)
	  
! Save output on this face
            face_output_field(output_face_number,1,1:6)=field(1:6)	  

! field on side 2, (xmax)
! E=-V/dl	    
	    field(Ex)=0d0
	    field(Ey)=-(Vy_positive_x2+Vy_negative_x2)/dl
	    field(Ez)=-(Vz_positive_x2+Vz_negative_x2)/dl
! H= I/dl	    
	    field(Hx)=0d0
	    field(Hy)= (Vz_positive_x2-Vz_negative_x2)/(Z0*dl)
	    field(Hz)=-(Vy_positive_x2-Vy_negative_x2)/(Z0*dl)
	  
! Save output on this face
            face_output_field(output_face_number,2,1:6)=field(1:6)	  
	    
	  end if		
	
	end if
