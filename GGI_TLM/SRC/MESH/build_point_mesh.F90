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
!SUBROUTINE build_point_mesh
!
! NAME
!     SUBROUTINE build_point_mesh
!
! DESCRIPTION
!     build_point_mesh:
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
SUBROUTINE build_point_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: point_number

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3
integer		:: triangle_number
integer		:: face_count
logical		:: set_mesh_flag

! START

  CALL write_line('CALLED: build_point_mesh',0,output_to_screen_flag)

  do point_number=1,n_points

     CALL point_to_cell(problem_points(point_number)%point,problem_points(point_number)%cell)
     CALL point_to_face(problem_points(point_number)%point,problem_points(point_number)%face)
    
  end do ! next point number

  CALL write_line('FINISHED: build_point_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE build_point_mesh
