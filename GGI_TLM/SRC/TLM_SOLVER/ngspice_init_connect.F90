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
!     Loop over mesh and link GGI_TLM voltage pulses to the ngspice solution
!
!     
! COMMENTS
!     The loop structure must be exactly the same as connect
!     
!     **** WE NEED TO ELIMINIATE UNUSED VARAIBLES ****
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
USE ngspice_F90

IMPLICIT NONE

! local variables
 
integer ( c_int ) :: istat
character*80 :: command_string

  integer cx,cy,cz
  integer face_number
  
  real*8 Vx,Vy,Vz
  
  real*8 Vx_min,Vy_min,Vz_min
  real*8 Vx_max,Vy_max,Vz_max
  
  real*8 Ix_min,Iy_min,Iz_min
  real*8 Ix_max,Iy_max,Iz_max
  
  integer	:: first_nz,last_nz
  integer	:: pol

! diode update variables  
  real*8	:: Is,nVt,Rs,Vd,Id,Zc,Vc,Ic
  integer	:: sign,diode_filter_number
  
! Spice circuit simulation link update vcariables
  real*8        :: Vspice
  integer       :: spice_node
   
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
  
  real*8	:: Is_min(3)
  real*8	:: Is_max(3)
  
  integer	:: hard_source_factor_min(6)
  integer	:: hard_source_factor_max(6)
  
  logical min_face_excitation
  
! Cable parameters  
    
  integer cable_face_junction_number
  
! START
  
  CALL write_line('CALLED: ngspice_init_connect',0,timestepping_output_to_screen_flag)
	
! Initial loop through the mesh getting the incident voltages for the ngspice solution

  face_number=0
  
! faces normal to x i.e. the y z plane
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
      
        if (cx.NE.1) then
            
          face_number=face_number+1
            
#include "ngspice_connect_normal_to_x.F90"            

        end if  ! cx.NE.1
  
! faces normal to y i.e. the x z plane
      
        if (cy.NE.1) then
 	
          face_number=face_number+1
            
#include "ngspice_connect_normal_to_y.F90"            
	            
        end if  ! cy.NE.1
      
        if (cz.NE.1) then
  
! faces normal to z i.e. the x y plane
      	
          face_number=face_number+1
            
#include "ngspice_connect_normal_to_z.F90"            
	
        end if
            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
          
  
  CALL write_line('FINISHED: ngspice_init_connect',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE ngspice_init_connect
