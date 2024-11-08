!
!    GGI_TLM Time domain electromagnetic field solver based on the TLM method
!    Copyright (C) 2013  Chris Smartt
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.   
!
! SUBROUTINE ngspice_init_connect
!
! NAME
!     ngspice_init_connect
!
! DESCRIPTION
!     Loop over mesh and work out the information to link GGI_TLM voltage pulses to the ngspice solution
!
!     
! COMMENTS
!     The loop structure must be exactly the same as connect
!     
!
! HISTORY
!
!     started  12/03/2019		CJS: Implement SPICE circuit model link
!
!
SUBROUTINE ngspice_init_connect

USE TLM_general
USE TLM_excitation
USE TLM_output
USE TLM_surface_materials
USE mesh
USE constants
USE file_information
USE iso_c_binding
#if defined(INCLUDE_NGSPICE)    
USE ngspice_F90
#endif

IMPLICIT NONE

! local variables
 
integer ( c_int ) :: istat
character*80 :: command_string

  integer cx,cy,cz
  integer face_number
    
! Spice circuit simulation link update vcariables
  integer       :: spice_port_count
  integer       :: spice_port_number
  integer       :: spice_node
   
! Material update parameters  
  integer material_type
  integer material_number
  
  integer,allocatable :: spice_port_set(:)
  
  integer :: i
  integer :: port,n1,n2
  
! START
  
  CALL write_line('CALLED: ngspice_init_connect',0,timestepping_output_to_screen_flag)
  
#if defined(INCLUDE_NGSPICE)    

! Initial loop through the mesh getting the incident voltages for the ngspice solution

  face_number=0
  spice_port_count=0
  
  ALLOCATE( spice_port_set(1:n_spice_ports) )
  spice_port_set(1:n_spice_ports)=0
  
! faces normal to x i.e. the y z plane
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
      
        if (cx.NE.1) then
            
         face_number=face_number+1
            
	
	 if (face_update_code(face_number).NE.0) then   ! not free space

!         face_update_code points to arrays which tell us about materials, excitations and outputs
	
	  material_number=abs( face_update_code_to_material_data(face_update_code(face_number),1) )
	  if (material_number.NE.0) then
	    material_type=surface_material_list(material_number)%type
	  else
	    material_type=0
	  end if
          	  						
	  if (material_type.NE.0) then ! not free space scatter
	
	  
            if (material_type.EQ.surface_material_type_SPICE) then	
                
              spice_port_number=surface_material_list(material_number)%Spice_circuit_file_port                
              if (spice_port_set(spice_port_number).NE.0) then
                write(*,*)'ERROR in ngspice_init_connect'
                write(*,*)'Spice port number',spice_port_number,' has already been set in the mesh'
                write(*,*)'Does the Spice surface cover more than one TLM cell face?'
                STOP 1 
              end if
             
	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-y').OR.	&
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+y') ) then
  
                spice_port_count=spice_port_count+1
                
                ng_material_number(spice_port_count)=material_number
                ng_sign(spice_port_count)=surface_material_list(material_number)%Spice_port_sign
                ng_spice_port(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_port
                ng_spice_node1(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
                ng_spice_node2(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(2)
                ng_face1(spice_port_count)=Vy_xmin
                ng_cx1(spice_port_count)=cx
                ng_cy1(spice_port_count)=cy
                ng_cz1(spice_port_count)=cz
                ng_face2(spice_port_count)=Vy_xmax
                ng_cx2(spice_port_count)=cx-1
                ng_cy2(spice_port_count)=cy
                ng_cz2(spice_port_count)=cz
                              	      	    
 	      else if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-z').OR.  &
	               (surface_material_list(material_number)%Spice_port_direction.EQ.'+z') ) then
  
                spice_port_count=spice_port_count+1
                
                ng_material_number(spice_port_count)=material_number
                ng_sign(spice_port_count)=surface_material_list(material_number)%Spice_port_sign
                ng_spice_port(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_port
                ng_spice_node1(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
                ng_spice_node2(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(2)
                ng_face1(spice_port_count)=Vz_xmin
                ng_cx1(spice_port_count)=cx
                ng_cy1(spice_port_count)=cy
                ng_cz1(spice_port_count)=cz
                ng_face2(spice_port_count)=Vz_xmin
                ng_cx2(spice_port_count)=cx-1
                ng_cy2(spice_port_count)=cy
                ng_cz2(spice_port_count)=cz
                      
              else
              
                write(*,*)'ERROR:'
                write(*,*)'Found SPICE face at cell ',cx,cy,cz,' face= face_xmin'
                write(*,*)'Spice port direction should be -y, +y, -z or +z'
                STOP 1
                             	      
              end if
              
              spice_port_set(spice_port_number)=1
           
	    end if  ! Spice surface type
	    
	  end if  ! not free space scatter

	end if  ! not free space

        end if  ! cx.NE.1
  
! faces normal to y i.e. the x z plane
      
        if (cy.NE.1) then
 	
          face_number=face_number+1
            
	
	if (face_update_code(face_number).NE.0) then   ! not free space

!         face_update_code points to arrays which tell us about materials, excitations and outputs
	
	  material_number=abs( face_update_code_to_material_data(face_update_code(face_number),1) )
	  if (material_number.NE.0) then
	    material_type=surface_material_list(material_number)%type
	  else
	    material_type=0
	  end if
          	  						
	  if (material_type.NE.0) then ! not free space scatter
	
	  
            if (material_type.EQ.surface_material_type_SPICE) then
                
              spice_port_number=surface_material_list(material_number)%Spice_circuit_file_port                
              if (spice_port_set(spice_port_number).NE.0) then
                write(*,*)'ERROR in ngspice_init_connect'
                write(*,*)'Spice port number',spice_port_number,' has already been set in the mesh'
                write(*,*)'Does the Spice surface cover more than one TLM cell face?'
                STOP 1 
              end if
            	 
	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-x').OR.	&
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+x') ) then
  
                spice_port_count=spice_port_count+1
                
                ng_material_number(spice_port_count)=material_number
                ng_sign(spice_port_count)=surface_material_list(material_number)%Spice_port_sign
                ng_spice_port(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_port
                ng_spice_node1(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
                ng_spice_node2(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(2)
                ng_face1(spice_port_count)=Vx_ymin
                ng_cx1(spice_port_count)=cx
                ng_cy1(spice_port_count)=cy
                ng_cz1(spice_port_count)=cz
                ng_face2(spice_port_count)=Vx_ymax
                ng_cx2(spice_port_count)=cx
                ng_cy2(spice_port_count)=cy-1
                ng_cz2(spice_port_count)=cz
                            	      
              else if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-z').OR.  &
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+z') ) then
  
                spice_port_count=spice_port_count+1
                
                ng_material_number(spice_port_count)=material_number
                ng_sign(spice_port_count)=surface_material_list(material_number)%Spice_port_sign
                ng_spice_port(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_port
                ng_spice_node1(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
                ng_spice_node2(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(2)
                ng_face1(spice_port_count)=Vz_ymin
                ng_cx1(spice_port_count)=cx
                ng_cy1(spice_port_count)=cy
                ng_cz1(spice_port_count)=cz
                ng_face2(spice_port_count)=Vz_ymax
                ng_cx2(spice_port_count)=cx
                ng_cy2(spice_port_count)=cy-1
                ng_cz2(spice_port_count)=cz
                      
              else
              
                write(*,*)'ERROR:'
                write(*,*)'Found SPICE face at cell ',cx,cy,cz,' face= face_ymin'
                write(*,*)'Spice port direction should be -x, +x, -z or +z'
                STOP 1
                             	      
                              	      
               end if
              
              spice_port_set(spice_port_number)=1
	                
	    end if  ! Spice surface type
	    
	  end if  ! not free space scatter

	end if  ! not free space
	            
        end if  ! cy.NE.1
      
        if (cz.NE.1) then
  
! faces normal to z i.e. the x y plane
      	
          face_number=face_number+1
            
	
	if (face_update_code(face_number).NE.0) then   ! not free space

!         face_update_code points to arrays which tell us about materials, excitations and outputs
	
	  material_number=abs( face_update_code_to_material_data(face_update_code(face_number),1) )
	  if (material_number.NE.0) then
	    material_type=surface_material_list(material_number)%type
	  else
	    material_type=0
	  end if
	  						
	  if (material_type.NE.0) then ! not free space scatter
		  
            if (material_type.EQ.surface_material_type_SPICE) then	 
                
              spice_port_number=surface_material_list(material_number)%Spice_circuit_file_port                
              if (spice_port_set(spice_port_number).NE.0) then
                write(*,*)'ERROR in ngspice_init_connect'
                write(*,*)'Spice port number',spice_port_number,' has already been set in the mesh'
                write(*,*)'Does the Spice surface cover more than one TLM cell face?'
                STOP 1 
              end if
              
	      if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-x').OR.	&
	          (surface_material_list(material_number)%Spice_port_direction.EQ.'+x') ) then
  
                spice_port_count=spice_port_count+1
                
                ng_material_number(spice_port_count)=material_number
                ng_sign(spice_port_count)=surface_material_list(material_number)%Spice_port_sign
                ng_spice_port(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_port
                ng_spice_node1(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
                ng_spice_node2(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(2)
                ng_face1(spice_port_count)=Vx_zmin
                ng_cx1(spice_port_count)=cx
                ng_cy1(spice_port_count)=cy
                ng_cz1(spice_port_count)=cz
                ng_face2(spice_port_count)=Vx_zmax
                ng_cx2(spice_port_count)=cx
                ng_cy2(spice_port_count)=cy
                ng_cz2(spice_port_count)=cz-1
                              	      
               else if ((surface_material_list(material_number)%Spice_port_direction.EQ.'-y').OR.	&
	           (surface_material_list(material_number)%Spice_port_direction.EQ.'+y') ) then
  
                spice_port_count=spice_port_count+1
                
                ng_material_number(spice_port_count)=material_number
                ng_sign(spice_port_count)=surface_material_list(material_number)%Spice_port_sign
                ng_spice_port(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_port
                ng_spice_node1(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(1)
                ng_spice_node2(spice_port_count)=surface_material_list(material_number)%Spice_circuit_file_nodes(2)
                ng_face1(spice_port_count)=Vy_zmin
                ng_cx1(spice_port_count)=cx
                ng_cy1(spice_port_count)=cy
                ng_cz1(spice_port_count)=cz
                ng_face2(spice_port_count)=Vy_zmax
                ng_cx2(spice_port_count)=cx
                ng_cy2(spice_port_count)=cy
                ng_cz2(spice_port_count)=cz-1
                      
              else
              
                write(*,*)'ERROR:'
                write(*,*)'Found SPICE face at cell ',cx,cy,cz,' face= face_zmin'
                write(*,*)'Spice port direction should be -x, +x, -y or +y'
                STOP 1                             	      

               end if	    
              
              spice_port_set(spice_port_number)=1
            
	     end if  ! Spice surface type
	    
	   end if  ! not free space scatter

	 end if  ! not free space
	
        end if
            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
  if (spice_port_count.NE.n_spice_ports) then
    write(*,*)'ERROR finding the spice ports in the mesh in ngspice_init_connect'
    write(*,*)'spice_port_count=',spice_port_count,' n_spice_ports=',n_spice_ports
    STOP 1
  end if
  
  do spice_port_count=1,n_spice_ports
  
    if (spice_port_set(spice_port_count).NE.1) then    
      write(*,*)'ERROR in ngspice_init_connect'
      write(*,*)'spice_port_count=',spice_port_count,' has not been found in the mesh'
      STOP 1     
    end if
    
    port=ng_spice_port(spice_port_count)
    n1=ng_spice_node1(spice_port_count)
    n2=ng_spice_node2(spice_port_count)
    
    write(*,*)'Checking port:'
    write(*,*)'From input file  : port=',port,' node1=',n1,' node2=',n2
    write(*,*)'From circuit file: port=',port,' node1=',cir_port_to_node1(port),' node2=',cir_port_to_node2(port)
    
    if ( (n1.EQ.cir_port_to_node1(port)).AND.(n2.EQ.cir_port_to_node2(port)) ) then
! we have the two nodes identified in the correct orientation
       
       write(*,*)'Port definition is consistent between input file and Spice_TEMPLATE.cir file'
       
    else
    
      write(*,*)'ERROR :Port definition is NOT consistent between input file and Spice_TEMPLATE.cir file'
      
      write(*,*)'WARNING :Port definition is NOT consistent between input file and Spice_TEMPLATE.cir file'
      write(*,*)'From input file  : port=',port,' node1=',n1,' node2=',n2
      write(*,*)'From circuit file: port=',port,' node1=',cir_port_to_node1(port),' node2=',cir_port_to_node2(port)
      
    end if
  
  end do
  
  DEALLOCATE( spice_port_set )
            
#endif
  
  CALL write_line('FINISHED: ngspice_init_connect',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE ngspice_init_connect
