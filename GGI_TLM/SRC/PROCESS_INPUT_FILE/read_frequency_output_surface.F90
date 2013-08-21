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
! SUBROUTINE read_frequency_output_surface
!
! NAME
!     read_frequency_output_surface
!
! DESCRIPTION
!     read frequency output surface packet
!
! Example packet:
!
!Frequency_output_surface_list
!1    		! number of Frequency output surfaces
!1		! FREQUENCY OUTPUT SURFACE NUMBER
!1     surface number
!1		! side of surface for output (+1 or -1)
!1.3E7 frequency for output
!Ex
!
! COMMENTS
!     
!
! HISTORY
!
!     started 5/12/2012 CJS
!
!
SUBROUTINE read_frequency_output_surface

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

! local variables

  integer	:: i
  integer	:: read_number
  integer	:: side_of_surface_for_output
  
! START  

  CALL write_line('CALLED: read_frequency_output_surface',0,output_to_screen_flag)
  
  read(input_file_unit,*,err=9000)n_frequency_output_surfaces
  
  CALL write_line_integer('number of frequency output surfaces',n_frequency_output_surfaces,0,output_to_screen_flag)
  
  if ( allocated( frequency_output_surface ) ) GOTO 9010
  
  ALLOCATE( frequency_output_surface(1:n_frequency_output_surfaces) )
  
  do i=1,n_frequency_output_surfaces
  
    CALL write_line_integer('Reading requency_output_surface number',i,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9000)read_number
    if (read_number.ne.i) GOTO 9020
  
    read(input_file_unit,*,err=9000)frequency_output_surface(i)%surface_number
    
    read(input_file_unit,*,err=9000)side_of_surface_for_output
    
    if (side_of_surface_for_output.eq.1) then
      frequency_output_surface(i)%output_on_outward_normal=.TRUE.
    else if (side_of_surface_for_output.eq.-1) then
      frequency_output_surface(i)%output_on_outward_normal=.FALSE.
    else 
      GOTO 9030
    end if
    
    read(input_file_unit,*,err=9000)frequency_output_surface(i)%frequency
      	
    CALL read_field_component_frequency_surface_output(input_file_unit,frequency_output_surface(i)%field_component)
     
  end do ! next frequency output surface    
  
  CALL write_line('FINISHED: read_frequency_output_surface',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading frequency_output_surface packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9010 CALL write_line('Error allocating frequency_output_surface:',0,.TRUE.)
     CALL write_line('frequency_output_surface already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading frequency_output_surface packet',0,.TRUE.)
     CALL write_line('frequency output surfaces should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9030 CALL write_line('Error reading frequency_output_surface packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_frequency_output_surface
