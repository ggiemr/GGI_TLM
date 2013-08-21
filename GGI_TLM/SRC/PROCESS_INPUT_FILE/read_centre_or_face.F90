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
! SUBROUTINE read_centre_or_face
!
! NAME
!     read_centre_or_face
!
! DESCRIPTION
!     read the point within a cell i.e. either the cell centre point or a face
!     this is used in specifying excitation and output points
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 13/09/2012 CJS
!
!
SUBROUTINE read_centre_or_face(file_unit,centre_or_face)

USE cell_parameters

IMPLICIT NONE

  integer		:: file_unit
  integer		:: centre_or_face

! local variables

character*256	:: input_line

! START  
    read(file_unit,'(A)',err=9000)input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,256)
    
    if      (input_line.eq.xmin_string) then   
      centre_or_face=face_xmin
    else if (input_line.eq.xmax_string) then   
      centre_or_face=face_xmax
    else if (input_line.eq.ymin_string) then   
      centre_or_face=face_ymin
    else if (input_line.eq.ymax_string) then   
      centre_or_face=face_ymax
    else if (input_line.eq.zmin_string) then   
      centre_or_face=face_zmin
    else if (input_line.eq.zmax_string) then   
      centre_or_face=face_zmax
    else if (input_line.eq.centre_string) then   
      centre_or_face=centre 
    else
    
      GOTO 9010
      
    end if

    
  RETURN
  
9000 CALL write_line('Error reading field component from file:',0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
     
9010 CALL write_line('Error reading field component',0,.TRUE.)
     CALL write_line("Expecting either 'centre', 'xmin', 'xmax', 'ymin', 'ymax' , 'zmin', 'zmax'",0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
  
END SUBROUTINE read_centre_or_face
