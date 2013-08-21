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
! SUBROUTINE connect
!
! NAME
!     connect
!
! DESCRIPTION
!     Loop over mesh and carry out the connection process.
!     faces normal to x are done first then those normal to y 
!     then those normal to z
!
!     All connections are applied on the xmin, ymin, zmin side of the cell,
!     A negative material number indicates that a material should be reveresed
!
!     Field outputs are saved on both sides of a face
!     Face output fields are calculated from incident and scattered voltage pulses
!     so as to be compatible with any materials or excitations present on the face
!     
!     The connection process applied is that approriate for free space,
!     PEC, PMC or thin layer material as appropriate
!
!
!     
! COMMENTS
!     Include a check for materials/ outputs at this stage,
!     Implementation of thin layer materials is included here and
!     Field output also becomes straightforward if recorded at this time
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!     Take account of surface normal and reverse surface material if required 11/09/2012 CJS
!     Add output on faces 13/09/2012 CJS
!     Start to include excitations on surfaces 18/09/2012 CJS
!     parallel 23/11/2012 CJS
!
!
SUBROUTINE connect

USE TLM_general
USE TLM_excitation
USE TLM_output
USE TLM_surface_materials
USE mesh
USE constants
USE file_information

IMPLICIT NONE

! local variables

  integer cx,cy,cz
  integer face_number
  
  real*8 Vx,Vy,Vz
  
  real*8 Vx_min,Vy_min,Vz_min
  real*8 Vx_max,Vy_max,Vz_max
  
  real*8 Ix_min,Iy_min,Iz_min
  real*8 Ix_max,Iy_max,Iz_max
  
  integer	:: first_nz,last_nz
   
! Material update parameters  
  integer material_type
  integer material_number
  logical reverse_material
  integer surface_filter_number
 
! Output parameters  
  
  integer output_number
  logical min_face_output
  
  integer output_face_number
  real*8 Vx_positive_y1,Vx_negative_y1
  real*8 Vx_positive_z1,Vx_negative_z1
  
  real*8 Vy_positive_x1,Vy_negative_x1
  real*8 Vy_positive_z1,Vy_negative_z1
  
  real*8 Vz_positive_x1,Vz_negative_x1
  real*8 Vz_positive_y1,Vz_negative_y1
  
  real*8 Vx_positive_y2,Vx_negative_y2
  real*8 Vx_positive_z2,Vx_negative_z2
  
  real*8 Vy_positive_x2,Vy_negative_x2
  real*8 Vy_positive_z2,Vy_negative_z2
  
  real*8 Vz_positive_x2,Vz_negative_x2
  real*8 Vz_positive_y2,Vz_negative_y2
  
  real*8	:: field(6)
  
! Excitation parameters  
    
  integer excitation_number
  integer excitation_face_number
  real*8	:: field_min(6)
  real*8	:: field_max(6)
  
  real*8	:: Js_min(3)
  real*8	:: Js_max(3)
  real*8	:: Ms_min(3)
  real*8	:: Ms_max(3)
  
  logical min_face_excitation
  
! Cable parameters  
    
  integer cable_face_junction_number

! START
  
  CALL write_line('CALLED: connect',0,timestepping_output_to_screen_flag)
	
! Check for any additional processes at each face and if there are:
! Add any thin layer models
! Save any requested outputs

  face_number=0
  excitation_face_number=0
  output_face_number=0
  
! faces normal to x i.e. the y z plane
  do cz=nz1,nz2
    do cy=1,ny
      do cx=2,nx
            
        face_number=face_number+1
	
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
	  
	    else if (material_type.EQ.surface_material_type_DISPERSIVE) then
	    
! Vy polarisation 	      
	      surface_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)
	      
              call surface_material_update(V(Vy_xmax,cx-1,cy,cz),Z0,V(Vy_xmin,cx,cy,cz),Z0,	&
	                                   Vy_max,Vy_min,material_number,reverse_material,surface_filter_number)	
	      
! Vz polarisation 	      
              call surface_material_update(V(Vz_xmax,cx-1,cy,cz),Z0,V(Vz_xmin,cx,cy,cz),Z0,	&
	                                   Vz_max,Vz_min,material_number,reverse_material,surface_filter_number+1)	
	    
	    end if	  
	    
	  end if	  

! Excitation		
          if (excitation_number.ne.0) then
	  
	    excitation_face_number=excitation_face_number+1
	    field_min(1:6)=face_excitation_field(excitation_face_number,1,1:6)
	    field_max(1:6)=face_excitation_field(excitation_face_number,2,1:6)
	    
	    Js_min(1:3)=face_excitation_field(excitation_face_number,1,7:9)
	    Js_max(1:3)=face_excitation_field(excitation_face_number,2,7:9)
	    Ms_min(1:3)=face_excitation_field(excitation_face_number,1,10:12)
	    Ms_max(1:3)=face_excitation_field(excitation_face_number,2,10:12)
	    
            Vy_min=Vy_min-field_min(Ey)*dl-field_min(Hz)*Z0*dl-Js_min(2)*dl*Z0/2d0-Ms_min(3)*dl/2d0
	    Vy_max=Vy_max-field_max(Ey)*dl+field_min(Hz)*Z0*dl-Js_max(2)*dl*Z0/2d0+Ms_max(3)*dl/2d0
	    
            Vz_min=Vz_min-field_min(Ez)*dl+field_min(Hy)*Z0*dl-Js_min(3)*dl*Z0/2d0+Ms_min(2)*dl/2d0
	    Vz_max=Vz_max-field_max(Ez)*dl-field_min(Hy)*Z0*dl-Js_min(3)*dl*Z0/2d0-Ms_min(2)*dl/2d0
	    
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
		
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
! faces normal to y i.e. the x z plane
  do cz=nz1,nz2
    do cy=2,ny
      do cx=1,nx
 	
        face_number=face_number+1
	
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
	  
	    else if (material_type.EQ.surface_material_type_DISPERSIVE) then
	    
	      surface_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)
	      
! Vx polarisation 	      
              call surface_material_update(V(Vx_ymax,cx,cy-1,cz),Z0,V(Vx_ymin,cx,cy,cz),Z0,	&
	                                   Vx_max,Vx_min,material_number,reverse_material,surface_filter_number)	
	      	      
! Vz polarisation 	      
              call surface_material_update(V(Vz_ymax,cx,cy-1,cz),Z0,V(Vz_ymin,cx,cy,cz),Z0,	&
	                                   Vz_max,Vz_min,material_number,reverse_material,surface_filter_number+1)	
					       	  
	    end if
	    
	  end if

! Excitation		
          if (excitation_number.ne.0) then
	  
	    excitation_face_number=excitation_face_number+1
	    field_min(1:6)=face_excitation_field(excitation_face_number,1,1:6)
	    field_max(1:6)=face_excitation_field(excitation_face_number,2,1:6)
	
	    Js_min(1:3)=face_excitation_field(excitation_face_number,1,7:9)
	    Js_max(1:3)=face_excitation_field(excitation_face_number,2,7:9)
	    Ms_min(1:3)=face_excitation_field(excitation_face_number,1,10:12)
	    Ms_max(1:3)=face_excitation_field(excitation_face_number,2,10:12)
	    
            Vx_min=Vx_min-field_min(Ex)*dl+field_min(Hz)*Z0*dl-Js_min(1)*dl*Z0/2d0+Ms_min(3)*dl/2d0
	    Vx_max=Vx_max-field_max(Ex)*dl-field_min(Hz)*Z0*dl-Js_max(1)*dl*Z0/2d0-Ms_max(3)*dl/2d0
	    
            Vz_min=Vz_min-field_min(Ez)*dl-field_min(Hx)*Z0*dl-Js_min(3)*dl*Z0/2d0-Ms_min(1)*dl/2d0
	    Vz_max=Vz_max-field_max(Ez)*dl+field_min(Hx)*Z0*dl-Js_max(3)*dl*Z0/2d0+Ms_max(1)*dl/2d0
	  
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
            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
! faces normal to z i.e. the x y plane

  if (rank.eq.0) then
    first_nz=2
  else
    first_nz=nz1
  end if
  
  last_nz=nz2

  do cz=first_nz,last_nz
    do cy=1,ny
      do cx=1,nx
      	
        face_number=face_number+1
	
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
	     
	    else if (material_type.EQ.surface_material_type_DISPERSIVE) then
	    
	      surface_filter_number=face_update_code_to_material_data(face_update_code(face_number),2)

! Vx polarisation 	      
              call surface_material_update(V(Vx_zmax,cx,cy,cz-1),Z0,V(Vx_zmin,cx,cy,cz),Z0,	&
	                                   Vx_max,Vx_min,material_number,reverse_material,surface_filter_number)	
	      	      
! Vy polarisation 	      
              call surface_material_update(V(Vy_zmax,cx,cy,cz-1),Z0,V(Vy_zmin,cx,cy,cz),Z0,	&
	                                   Vy_max,Vy_min,material_number,reverse_material,surface_filter_number+1)	
	  
	    end if
	    
	  end if

! Excitation		
          if (excitation_number.ne.0) then
	  
	    excitation_face_number=excitation_face_number+1
	    field_min(1:6)=face_excitation_field(excitation_face_number,1,1:6)
	    field_max(1:6)=face_excitation_field(excitation_face_number,2,1:6)

	    Js_min(1:3)=face_excitation_field(excitation_face_number,1,7:9)
	    Js_max(1:3)=face_excitation_field(excitation_face_number,2,7:9)
	    Ms_min(1:3)=face_excitation_field(excitation_face_number,1,10:12)
	    Ms_max(1:3)=face_excitation_field(excitation_face_number,2,10:12)
	    
            Vx_min=Vx_min-field_min(Ex)*dl-field_min(Hy)*Z0*dl-Js_min(1)*dl*Z0/2d0-Ms_min(2)*dl/2d0
	    Vx_max=Vx_max-field_max(Ex)*dl+field_min(Hy)*Z0*dl-Js_max(1)*dl*Z0/2d0+Ms_max(2)*dl/2d0
	    
            Vy_min=Vy_min-field_min(Ey)*dl+field_min(Hx)*Z0*dl-Js_min(2)*dl*Z0/2d0+Ms_min(1)*dl/2d0
	    Vy_max=Vy_max-field_max(Ey)*dl-field_min(Hx)*Z0*dl-Js_max(2)*dl*Z0/2d0-Ms_max(1)*dl/2d0
	  
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
            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
  CALL write_line('FINISHED: connect',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE connect
