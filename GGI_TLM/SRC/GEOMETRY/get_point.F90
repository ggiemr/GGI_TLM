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
!SUBROUTINE get_point
!
! NAME
!     SUBROUTINE get_point
!
! DESCRIPTION
!     get_point:
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE get_point(point_number,cell,point)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

  integer		:: point_number
  type(ijk) 		:: cell
  type(xyz)		:: point

! local variables

! START

  if ( (point_number.lt.1).OR.(point_number.gt.n_points) ) GOTO 9000

  cell%i=problem_points(point_number)%cell%i
  cell%j=problem_points(point_number)%cell%j
  cell%k=problem_points(point_number)%cell%k
  
  point%x=problem_points(point_number)%point%x
  point%y=problem_points(point_number)%point%y
  point%z=problem_points(point_number)%point%z
  
  RETURN
  
9000 CALL write_line('Error in get_point:',0,.TRUE.)
     CALL write_line('Point number does not exist',0,.TRUE.)
     CALL write_line_integer('Point number requested',point_number,0,.TRUE.)
     CALL write_line_integer('Number of points',n_points,0,.TRUE.)
     STOP

  
END SUBROUTINE get_point
