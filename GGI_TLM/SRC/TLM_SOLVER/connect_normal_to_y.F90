	
	if (face_update_code(face_number).eq.0) then   ! free space update, no excitation, no output
	
	
          Vx=V(Vx_ymin,cx,cy,cz)+ V(Vx_ymax,cx,cy-1,cz)
	  V(Vx_ymin,cx,cy,cz)  =  Vx-V(Vx_ymin,cx,cy,cz)
          V(Vx_ymax,cx,cy-1,cz)=  Vx-V(Vx_ymax,cx,cy-1,cz)
	
          Vz=V(Vz_ymin,cx,cy,cz)+ V(Vz_ymax,cx,cy-1,cz)
	  V(Vz_ymin,cx,cy,cz)  =  Vz-V(Vz_ymin,cx,cy,cz)
          V(Vz_ymax,cx,cy-1,cz)=  Vz-V(Vz_ymax,cx,cy-1,cz)
	  
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
	    Vx_negative_y1=V(Vx_ymin,cx,cy,cz)
	    Vz_negative_y1=V(Vz_ymin,cx,cy,cz)
	    Vx_positive_y2=V(Vx_ymax,cx,cy-1,cz)
	    Vz_positive_y2=V(Vz_ymax,cx,cy-1,cz)
	    
	  end if		
		
	  if (material_type.EQ.0) then ! no material, simple connect
	
            Vx_min=V(Vx_ymin,cx,cy,cz)+ V(Vx_ymax,cx,cy-1,cz)
	    Vx_max=Vx_min
	
            Vz_min=V(Vz_ymin,cx,cy,cz)+ V(Vz_ymax,cx,cy-1,cz)
	    Vz_max=Vz_min
	  
	  else ! material_type.NE.0
	  
	    if (material_type.EQ.surface_material_type_PEC) then
	  
              Vx_min=0D0
	      Vx_max=0D0
	
              Vz_min=0D0
	      Vz_max=0D0
	  
	    else if (material_type.EQ.surface_material_type_PMC) then
	  
              Vx_min=2D0*V(Vx_ymin,cx,cy,cz)
	      Vx_max=2D0*V(Vx_ymax,cx,cy-1,cz)
	
              Vz_min=2D0*V(Vz_ymin,cx,cy,cz)
	      Vz_max=2D0*V(Vz_ymax,cx,cy-1,cz)
	  
	    else if (material_type.EQ.surface_material_type_FREE_SPACE) then	 
	     
              Vx_min=V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy-1,cz)
	      Vx_max=Vx_min
	
              Vz_min=V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy-1,cz)
	      Vz_max=Vz_min
	  
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
  	      
	        CALL diode_calc(V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy-1,cz),Z0/2d0,Is,nVt,Rs,Zc,Vc,Vd,Id,Ic,sign)
		Vx_min=Vd
		Vx_max=Vd
		
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       Ic)
	      
	      else ! diode is not is this direction so do free space update for x polarisation
	      
                Vx_min=V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy-1,cz)
	        Vx_max=Vx_min
       
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
	      
	        CALL diode_calc(V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy-1,cz),Z0/2d0,Is,nVt,Rs,Zc,Vc,Vd,Id,Ic,sign)
		Vz_min=Vd
		Vz_max=Vd
		
                CALL evaluate_Zfilter (surface_material_list(material_number)%Diode_Cj_Z,	 &
	        		       Diode_Cj_filter_data(diode_filter_number),	 &
			  	       Ic)
	      
	      else ! diode is not is this direction so do free space update for z polarisation
	      
                Vz_min=V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy-1,cz)
	        Vz_max=Vz_min
       
              end if
	  
	    else if ( (material_type.EQ.surface_material_type_DISPERSIVE).OR.			&
	              (material_type.EQ.surface_material_type_ANISOTROPIC_DISPERSIVE) ) then
	    
	      surface_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)
	      
! Vx polarisation 	      
              pol=1
              call surface_material_update(V(Vx_ymax,cx,cy-1,cz),Z0,V(Vx_ymin,cx,cy,cz),Z0,	&
	                                   Vx_max,Vx_min,material_number,pol,reverse_material,surface_filter_number)	
	      	      
! Vz polarisation 	      
              pol=3
              call surface_material_update(V(Vz_ymax,cx,cy-1,cz),Z0,V(Vz_ymin,cx,cy,cz),Z0,	&
	                                   Vz_max,Vz_min,material_number,pol,reverse_material,surface_filter_number+1)	
	  
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
	    
            Vx_min=Vx_min*hard_source_factor_min(Ex)-field_min(Ex)*dl+field_min(Hz)*Z0*dl+Js_min(1)*dl*Z0/2d0-Ms_min(3)*dl/2d0
	    Vx_max=Vx_max*hard_source_factor_max(Ex)-field_max(Ex)*dl-field_max(Hz)*Z0*dl+Js_max(1)*dl*Z0/2d0+Ms_max(3)*dl/2d0
	    
            Vz_min=Vz_min*hard_source_factor_min(Ez)-field_min(Ez)*dl-field_min(Hx)*Z0*dl+Js_min(3)*dl*Z0/2d0+Ms_min(1)*dl/2d0
	    Vz_max=Vz_max*hard_source_factor_max(Ez)-field_max(Ez)*dl+field_max(Hx)*Z0*dl+Js_max(3)*dl*Z0/2d0-Ms_max(1)*dl/2d0
	  
	  end if
	  					       
! complete the connection					       
          V(Vx_ymax,cx,cy-1,cz)  =  Vx_max-V(Vx_ymax,cx,cy-1,cz)
	  V(Vx_ymin,cx,cy,cz)	 =  Vx_min-V(Vx_ymin,cx,cy,cz)
	    				   
          V(Vz_ymax,cx,cy-1,cz)  =  Vz_max-V(Vz_ymax,cx,cy-1,cz)
	  V(Vz_ymin,cx,cy,cz)	 =  Vz_min-V(Vz_ymin,cx,cy,cz)

! Cable connect		
          if (cable_face_junction_number.ne.0) then
	    
	    call Face_Cable_Junction( cable_face_junction_number )
	  
	  end if
	    
! Output		
          if (output_number.ne.0) then
	  	  
! get the scattered voltage pulses required for output field calculation	  
	    Vx_positive_y1=V(Vx_ymin,cx,cy,cz)
	    Vz_positive_y1=V(Vz_ymin,cx,cy,cz)
	    Vx_negative_y2=V(Vx_ymax,cx,cy-1,cz)
	    Vz_negative_y2=V(Vz_ymax,cx,cy-1,cz)
	    
	    output_face_number=output_face_number+1
	    
! field on side 1 (ymin)

! E=-V/dl	    
	    field(Ex)=-(Vx_positive_y1+Vx_negative_y1)/dl
	    field(Ey)=0d0
	    field(Ez)=-(Vz_positive_y1+Vz_negative_y1)/dl
! H= I/dl	    
	    field(Hx)=-(Vz_positive_y1-Vz_negative_y1)/(Z0*dl)
	    field(Hy)=0d0
	    field(Hz)= (Vx_positive_y1-Vx_negative_y1)/(Z0*dl)
	  
! Save output on this face
            face_output_field(output_face_number,1,1:6)=field(1:6)	  
	    
! field on side 2 (ymax)

! E=-V/dl	    
	    field(Ex)=-(Vx_positive_y2+Vx_negative_y2)/dl
	    field(Ey)=0d0
	    field(Ez)=-(Vz_positive_y2+Vz_negative_y2)/dl
! H= I/dl	    
	    field(Hx)=-(Vz_positive_y2-Vz_negative_y2)/(Z0*dl)
	    field(Hy)=0d0
	    field(Hz)= (Vx_positive_y2-Vx_negative_y2)/(Z0*dl)
	  
! Save output on this face
            face_output_field(output_face_number,2,1:6)=field(1:6)	  
	    	  
	  end if ! output
	
	end if ! special cell_face
