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
! SUBROUTINE read_cable_output_list
!
! NAME
!     read_cable_output_list
!
! DESCRIPTION
!     read cable output list packet
!
! Example packet:
!
!Cable_output_list
!1  number of cable_outputs
!1  CABLE OUTPUT NUMBER
!1  cable number
!1  output point number (closest point)
!
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!
!
SUBROUTINE read_cable_output_list

USE TLM_general
USE Cables
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: cable_output_number
integer :: read_number

integer :: i

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_cable_output_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_cable_outputs
  
  CALL write_line_integer('number of cable outputs',n_cable_outputs,0,output_to_screen_flag)
  
  if ( allocated( cable_output_list ) ) GOTO 9000
  
  ALLOCATE ( cable_output_list(1:n_cable_outputs) )

  do cable_output_number=1,n_cable_outputs
  
    CALL write_line_integer('Reading cable number',cable_output_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.cable_output_number) goto 9010
 
! read the cable number
    read(input_file_unit,*,err=9005)cable_output_list(cable_output_number)%cable_number
 
! read the closest output point number
    read(input_file_unit,*,err=9005)cable_output_list(cable_output_number)%closest_point_number
      
  end do ! next cable output

  CALL write_line('FINISHED: read_cable_output_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating cable_output_list:',0,.TRUE.)
     CALL write_line('cable_output_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading cable_output_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading cable_output_list packet data',0,.TRUE.)
     CALL write_line('Cables should be numbered in order',0,.TRUE.)
     STOP

9020 CALL write_line('Error reading cable_output_list packet data',0,.TRUE.)
     CALL write_line_integer('Error reading the line list, cable=',cable_output_number,0,.TRUE.)
     STOP
  
  
END SUBROUTINE read_cable_output_list
