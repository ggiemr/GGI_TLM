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
!Far_field_surface_list
!2     number of far field surfaces
!1     FAR FIELD SURFACE NUMBER
!4     surface number
!1		! side of surface for output (+1 or -1)
!300E6 frequency for far field calculation
!0.0 180.0 1.0  Theta_min  Theta_max Theta_step
!0.0 90.0 90.0    Phi_min  Phi_max Phi_step
!2     FAR FIELD SURFACE NUMBER
!4     surface number
!1		! side of surface for output (+1 or -1)
!400E6 frequency for far field calculation
!0.0 180.0 1.0  Theta_min  Theta_max Theta_step
!0.0 90.0 90.0    Phi_min  Phi_max Phi_step
!
!Optional request for 2D far field transform includes 2D and a line with the 'z' axis direction
!2D
!z
!
! COMMENTS
!     
!
! HISTORY
!
!       started 5/12/2012 CJS
!	allow multiple far field surfaces 23/5/2014
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
  integer 	:: read_number
  integer 	:: surface
  character	:: ch
  character*2	:: ch2

! START  

  CALL write_line('CALLED: read_far_field_surface',0,output_to_screen_flag)
  
  if (n_far_field_surfaces.NE.0) GOTO 9000
  
  read(input_file_unit,*,err=9005,end=9005)n_far_field_surfaces
  
  ALLOCATE( far_field_surface(1:n_far_field_surfaces) )
  
  do surface=1,n_far_field_surfaces
  
    CALL write_line_integer('Reading far field surface number',surface,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005,end=9005)read_number
    if (read_number.ne.surface) goto 9007
  
    read(input_file_unit,*,err=9005,end=9005)far_field_surface(surface)%surface_number
    
    read(input_file_unit,*,err=9005,end=9005)side_of_surface_for_output
    
    if (side_of_surface_for_output.eq.1) then
      far_field_surface(surface)%output_on_outward_normal=.TRUE.
    else if (side_of_surface_for_output.eq.-1) then
      far_field_surface(surface)%output_on_outward_normal=.FALSE.
    else 
      GOTO 9010
    end if
    
    read(input_file_unit,*,err=9005,end=9005)far_field_surface(surface)%frequency
      
    read(input_file_unit,*,err=9005,end=9005)far_field_surface(surface)%theta_min,        &
  				    far_field_surface(surface)%theta_max,        &
  				    far_field_surface(surface)%theta_step
      
    read(input_file_unit,*,err=9005,end=9005)far_field_surface(surface)%phi_min,  &
  				    far_field_surface(surface)%phi_max,  &
                                    far_field_surface(surface)%phi_step
    
! convert angles to radians
    far_field_surface(surface)%theta_min =(pi/180d0)*far_field_surface(surface)%theta_min 
    far_field_surface(surface)%theta_max =(pi/180d0)*far_field_surface(surface)%theta_max 
    far_field_surface(surface)%theta_step=(pi/180d0)*far_field_surface(surface)%theta_step

    far_field_surface(surface)%phi_min =(pi/180d0)*far_field_surface(surface)%phi_min  
    far_field_surface(surface)%phi_max =(pi/180d0)*far_field_surface(surface)%phi_max  
    far_field_surface(surface)%phi_step=(pi/180d0)*far_field_surface(surface)%phi_step  

! assume 3D transformation is used unless specified otherwise  
    far_field_surface(surface)%dim=3
    far_field_surface(surface)%direction=3
  
! check for 2D transformation request
    
    read(input_file_unit,'(A2)',err=1000)ch2
    
    if ( (ch2.eq.'2D').OR.(ch2.eq.'2d') ) then
    
      far_field_surface(surface)%dim=2
      
      read(input_file_unit,'(A)',err=9020,end=9020)ch
      
      if ( (ch.eq.'x').OR.(ch.eq.'X') ) then
        far_field_surface(surface)%direction=1
      else if ( (ch.eq.'y').OR.(ch.eq.'Y') ) then
        far_field_surface(surface)%direction=2
      else if ( (ch.eq.'z').OR.(ch.eq.'Z') ) then
        far_field_surface(surface)%direction=3
      else 
        goto 9020
      end if
      
      if (far_field_surface(surface)%direction.NE.3) then
        write(*,*)'ERROR in read_far_field_surface'
	write(*,*)'Cannot have the 2d axis direction anything other than z at the moment'
	STOP
      end if

! for the 2D transform  we only step in the phi direction
      far_field_surface(surface)%theta_min =pi/2d0
      far_field_surface(surface)%theta_max =pi/2d0
      far_field_surface(surface)%theta_step=1d0
        
    else
! no 2D transform has been requested so step back in the input file. 
  
      backspace(unit=input_file_unit)
  
    end if
    
1000 CONTINUE

  end do ! next surface 
  
  CALL write_line('FINISHED: read_far_field_surface',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error: far_field_surface_list data already set',0,.TRUE.)
     STOP 1
    
9005 CALL write_line('Error reading far_field_surface_list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9007 CALL write_line('Error reading far_field_surface_list list packet',0,.TRUE.)
     CALL write_line('Far field surfaces should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
        
9010 CALL write_line('Error reading far_field_surface_list packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
       
9020 CALL write_line('Error reading far_field_surface_list packet',0,.TRUE.)
     CALL write_line("Problem with 2D far field transformation specification",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1

END SUBROUTINE read_far_field_surface
