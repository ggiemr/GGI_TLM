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
!SUBROUTINE build_line_geometry
!
! NAME
!     SUBROUTINE build_line_geometry
!
! DESCRIPTION
!     build_line_geometry:
!
!     Create line geometric entities from the defined line type
!     The geometric entities consist of line segments
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
SUBROUTINE build_line_geometry

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: line_number

! START

  CALL write_line('CALLED: Build_line_geometry',0,output_to_screen_flag)

  do line_number=1,n_lines

    if (problem_lines(line_number)%line_type.EQ.line_type_straight_line) then
      
      CALL build_line_straight_line(line_number)
      
    else if (problem_lines(line_number)%line_type.EQ.line_type_straight_line2) then
      
      CALL build_line_straight_line2(line_number)
       
    else if (problem_lines(line_number)%line_type.EQ.line_type_arc) then
      
      CALL build_line_arc(line_number)
       
    else ! line type not yet defined
    
      GOTO 9000
      
    end if ! line_type

  end do ! next line number
  
  CALL plot_line_segments()

  CALL write_line('FINISHED: Build_line_geometry',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error in build_line_geometry:',0,.TRUE.)
     CALL write_line('line type number not defined',0,.TRUE.)
     CALL write_line_integer('line type number',problem_lines(line_number)%line_type,0,.TRUE.)
     STOP

  
END SUBROUTINE build_line_geometry
