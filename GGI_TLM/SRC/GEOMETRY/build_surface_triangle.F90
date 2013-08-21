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
!SUBROUTINE build_surface_triangle
!
! NAME
!     SUBROUTINE build_surface_triangle
!
! DESCRIPTION
!     build_surface_triangle:
!
!     create a triangulated triangle surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 12/09/2012 CJS
!
!
SUBROUTINE build_surface_triangle(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3

! START

! CREATE A TRIANGULATED RECTANGULAR BLOCK DEFINED BY OPPOSITE CORNER COORDINATES
    
      number_of_triangles=1
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! create triangle vertices

      point1%x=problem_surfaces(surface_number)%surface_parameters(1)
      point1%y=problem_surfaces(surface_number)%surface_parameters(2)
      point1%z=problem_surfaces(surface_number)%surface_parameters(3)
      
      point2%x=problem_surfaces(surface_number)%surface_parameters(4)
      point2%y=problem_surfaces(surface_number)%surface_parameters(5)
      point2%z=problem_surfaces(surface_number)%surface_parameters(6)
      
      point3%x=problem_surfaces(surface_number)%surface_parameters(7)
      point3%y=problem_surfaces(surface_number)%surface_parameters(8)
      point3%z=problem_surfaces(surface_number)%surface_parameters(9)
    
! apply the transformation to each of the points
      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point2,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
	  
! create surface triangle
      triangle_count=0
	  
      triangle_count=triangle_count+1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point2
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3

  RETURN

END SUBROUTINE build_surface_triangle
