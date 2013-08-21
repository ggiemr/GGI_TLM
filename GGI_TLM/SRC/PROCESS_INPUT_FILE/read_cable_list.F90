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
! SUBROUTINE read_cable_list
!
! NAME
!     read_cable_list
!
! DESCRIPTION
!     read cable list packet
!
! Example packet:
!
!Cable_list
!2  number of cables
!1  	CABLE NUMBER
!1 	cable geometry number
!1 	number of lines on cable route
!1 	cable line list
!1  	end 1 junction number
!2  	end 2 junction number
!2  	CABLE NUMBER
!1 	cable geometry number
!3 	number of lines on cable route
!2 3 4	cable line list
!2  	end 2 junction number
!3  	end 3 junction number
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!
!
SUBROUTINE read_cable_list

USE TLM_general
USE Cables
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: cable_number
integer :: read_number

integer :: i

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_cable_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_cables
  
  CALL write_line_integer('number of cables',n_cables,0,output_to_screen_flag)
  
  if ( allocated( cable_list ) ) GOTO 9000
  
  ALLOCATE ( cable_list(1:n_cables) )

  do cable_number=1,n_cables
  
    CALL write_line_integer('Reading cable number',cable_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.cable_number) goto 9010
 
! read the cable geometry number 
    read(input_file_unit,*,err=9005)cable_list(cable_number)%cable_geometry_number
 
! read the number of lines on the cable route   
    read(input_file_unit,*,err=9005)cable_list(cable_number)%n_lines
    
    if (cable_list(cable_number)%n_lines.gt.0) then
! allocate and read the line list
    
      ALLOCATE ( cable_list(cable_number)%line_list(1:cable_list(cable_number)%n_lines) )

      read(input_file_unit,*,err=9020)(cable_list(cable_number)%line_list(i),i=1,cable_list(cable_number)%n_lines)
      
    end if   ! n_lines.GT.0
    
    read(input_file_unit,*,err=9005)cable_list(cable_number)%junction_1
    read(input_file_unit,*,err=9005)cable_list(cable_number)%junction_2
      
  end do ! next cable 

  CALL write_line('FINISHED: read_cable_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating cable_list:',0,.TRUE.)
     CALL write_line('cable_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading cable_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading cable_list packet data',0,.TRUE.)
     CALL write_line('Cables should be numbered in order',0,.TRUE.)
     STOP

9020 CALL write_line('Error reading cable_list packet data',0,.TRUE.)
     CALL write_line_integer('Error reading the line list, cable=',cable_number,0,.TRUE.)
     STOP
  
  
END SUBROUTINE read_cable_list
