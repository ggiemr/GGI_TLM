	
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
	    
            
	    end if  ! Spice link
	    
	  end if  ! not free space scatter

	end if  ! not free space
