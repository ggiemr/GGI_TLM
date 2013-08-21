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
! SUBROUTINE read_far_field_surface
!
! NAME
!     read_far_field_surface
!
! DESCRIPTION
!     read far field surface packet
!
! Example packet:
!
!Far_field_surface
!1     surface number
!1		! side of surface for output (+1 or -1)
!1.3E7 frequency for far field calculation
!0.0 180.0 1.0  Theta_min  Theta_max Theta_step
!0.0 90.0 90.0    Phi_min  Phi_max Phi_step
!
! COMMENTS
!     
!
! HISTORY
!
!     started 5/12/2012 CJS
!
!
SUBROUTINE read_far_field_surface

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

! local variables

  integer	:: side_of_surface_for_output

! START  

  CALL write_line('CALLED: read_far_field_surface',0,output_to_screen_flag)
  
  n_far_field_surfaces=1
  
  read(input_file_unit,*,err=9000)far_field_surface%surface_number
    
  read(input_file_unit,*,err=9000)side_of_surface_for_output
    
  if (side_of_surface_for_output.eq.1) then
    far_field_surface%output_on_outward_normal=.TRUE.
  else if (side_of_surface_for_output.eq.-1) then
    far_field_surface%output_on_outward_normal=.FALSE.
  else 
    GOTO 9010
  end if
    
  read(input_file_unit,*,err=9000)far_field_surface%frequency
      
  read(input_file_unit,*,err=9000)far_field_surface%theta_min,        &
  				  far_field_surface%theta_max,        &
  				  far_field_surface%theta_step
      
  read(input_file_unit,*,err=9000)far_field_surface%phi_min,  &
  				  far_field_surface%phi_max,  &
                                    far_field_surface%phi_step
    
! convert angles to radians
  far_field_surface%theta_min =(pi/180d0)*far_field_surface%theta_min 
  far_field_surface%theta_max =(pi/180d0)*far_field_surface%theta_max 
  far_field_surface%theta_step=(pi/180d0)*far_field_surface%theta_step

  far_field_surface%phi_min =(pi/180d0)*far_field_surface%phi_min  
  far_field_surface%phi_max =(pi/180d0)*far_field_surface%phi_max  
  far_field_surface%phi_step=(pi/180d0)*far_field_surface%phi_step  
  
  CALL write_line('FINISHED: read_far_field_surface',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading far_field_surface packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
       
9010 CALL write_line('Error reading far_field_surface packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

END SUBROUTINE read_far_field_surface
