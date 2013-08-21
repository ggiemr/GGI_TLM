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
! SUBROUTINE read_rcs_surface
!
! NAME
!     read_rcs_surface
!
! DESCRIPTION
!     read RCS surface packet
!
! Example packet:
!
!RCS_surface
!2     surface number
!1     side of surface for output
!0.1e8 3e8 0.01e8 frequency range for RCS calculation
!180.0  Theta 
!0.0  Phi
!
! COMMENTS
!     
!
! HISTORY
!
!     started 2/1/2013 CJS based on FieldSolve code
!
!
SUBROUTINE read_rcs_surface

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

! local variables

  integer	:: side_of_surface_for_output

! START  

  CALL write_line('CALLED: read_rcs_surface',0,output_to_screen_flag)
  
  n_rcs_surfaces=1
  
  read(input_file_unit,*,err=9000)rcs_surface%surface_number
    
  read(input_file_unit,*,err=9000)side_of_surface_for_output
    
  if (side_of_surface_for_output.eq.1) then
    rcs_surface%output_on_outward_normal=.TRUE.
  else if (side_of_surface_for_output.eq.-1) then
    rcs_surface%output_on_outward_normal=.FALSE.
  else 
    GOTO 9010
  end if
  
  read(input_file_unit,*,err=9000)rcs_surface%fmin,        &
  				  rcs_surface%fmax,        &
  				  rcs_surface%fstep
      
  read(input_file_unit,*,err=9000)rcs_surface%theta
      
  read(input_file_unit,*,err=9000)rcs_surface%phi
    
! convert angles to radians
  rcs_surface%theta =(pi/180d0)*rcs_surface%theta 

  rcs_surface%phi =(pi/180d0)*rcs_surface%phi  
  
  CALL write_line('FINISHED: read_rcs_surface',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading rcs_surface packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
       
9010 CALL write_line('Error reading rcs_surface packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

END SUBROUTINE read_rcs_surface
