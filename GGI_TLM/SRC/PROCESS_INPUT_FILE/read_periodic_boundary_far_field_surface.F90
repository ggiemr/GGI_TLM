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
! SUBROUTINE read_periodic_boundary_far_field_surface
!
! NAME
!     read_periodic_boundary_far_field_surface
!
! DESCRIPTION
!     read far field surface packet
!
! Example packet:
!
!Periodic_boundary_far_field_surface
!2 		! number of periodic far field surfaces
!1 		! PERIODIC FAR FIELD SURFACE NUMBER
!3 		! surface number
!1		! side of surface for output (1 or -1)
!0.1e9 18e9 10e6 ! fmin fmax fstep
!transmission
!0 0		! order of grating lobe in x and y
!2 		! PERIODIC FAR FIELD SURFACE NUMBER
!2 		! surface number
!1		! side of surface for output (1 or -1)
!1e9 18e9 10e6 	! fmin fmax fstep
!reflection
!0 0		! order of grating lobe in x and y
!
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 1/5/2014 CJS based on HIRF-SE code Fieldsolve
!
!
SUBROUTINE read_periodic_boundary_far_field_surface

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

! local variables

  integer	:: surface
  integer	:: PB_surface_number
  integer	:: side_of_surface_for_output
  character	:: ch

! START  

  CALL write_line('CALLED: read_periodic_boundary_far_field_surface',0,output_to_screen_flag)
  
  read(input_file_unit,*,err=9000)n_PB_far_field_surfaces
  
  CALL write_line_integer('number of periodic boundary far field surfaces',n_PB_far_field_surfaces,0,output_to_screen_flag)
 
  if ( allocated( PB_far_field_surface ) ) GOTO 9010
  
  ALLOCATE( PB_far_field_surface(1:n_PB_far_field_surfaces) )
  
  do surface=1,n_PB_far_field_surfaces
  
    read(input_file_unit,*,err=9000)PB_surface_number
    if (PB_surface_number.ne.surface) goto 9020
  
    read(input_file_unit,*,err=9000)PB_far_field_surface(surface)%surface_number
    
    read(input_file_unit,*,err=9000)side_of_surface_for_output
    
    if (side_of_surface_for_output.eq.1) then
      PB_far_field_surface(surface)%output_on_outward_normal=.TRUE.
    else if (side_of_surface_for_output.eq.-1) then
      PB_far_field_surface(surface)%output_on_outward_normal=.FALSE.
    else 
      GOTO 9030
    end if
    
    read(input_file_unit,*,err=9000)PB_far_field_surface(surface)%fmin,        &
  				    PB_far_field_surface(surface)%fmax,        &
  				    PB_far_field_surface(surface)%fstep
      
 ! check for reflection or transmission side of surface
    
    read(input_file_unit,'(A)',err=9000)ch
    
    if ( (ch.eq.'r').OR.(ch.eq.'R') ) then    
      PB_far_field_surface(surface)%r_t_option='r'     
    else if ( (ch.eq.'t').OR.(ch.eq.'T') ) then    
      PB_far_field_surface(surface)%r_t_option='t'     
    else     
      GOTO 9040  
    end if
    
    read(input_file_unit,*,err=9000)PB_far_field_surface(surface)%m,PB_far_field_surface(surface)%n
    
  end do ! next periodic_boundary_far_field_surfac
  
  CALL write_line('FINISHED: read_periodic_boundary_far_field_surface',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading periodic_boundary_far_field_surface packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
    
9010 CALL write_line('Error in periodic_boundary_far_field_surface packet data from input file:',0,.TRUE.)
     CALL write_line('far_field_surface packet data is already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error in periodic_boundary_far_field_surface packet data from input file:',0,.TRUE.)
     CALL write_line('Periodic boundary far field surfaces should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
       
9030 CALL write_line('Error reading periodic_boundary_far_field_surface packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
       
9040 CALL write_line('Error reading periodic_boundary_far_field_surface packet',0,.TRUE.)
     CALL write_line("Expecting either Reflection or Transmission to be specified...",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

END SUBROUTINE read_periodic_boundary_far_field_surface
