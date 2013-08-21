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
!SUBROUTINE build_surface_xplane
!
! NAME
!     SUBROUTINE build_surface_xplane
!
! DESCRIPTION
!     build_surface_xplane:
!
!     create a triangulated xplane surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_surface_xplane(surface_number)

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

real*8	:: x1,y1,z1,x2,y2,z2
real*8	:: swap

type(xyz)	:: point1,point2,point3,point4

! START

! CREATE A TRIANGULATED RECTANGULAR BLOCK DEFINED BY OPPOSITE CORNER COORDINATES
    
      number_of_triangles=2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! create rectangular_block vertices
      x1=problem_surfaces(surface_number)%surface_parameters(1)
      y1=problem_surfaces(surface_number)%surface_parameters(2)
      z1=problem_surfaces(surface_number)%surface_parameters(3)
      x2=problem_surfaces(surface_number)%surface_parameters(4)
      y2=problem_surfaces(surface_number)%surface_parameters(5)
      z2=problem_surfaces(surface_number)%surface_parameters(6)
      
! ensure that x1=x2

      if (x1.NE.x2) then
        GOTO 9000
      end if

      point1%x=x1
      point1%y=y1
      point1%z=z1
      point2%x=x1
      point2%y=y2
      point2%z=z1
      point3%x=x1
      point3%y=y2
      point3%z=z2
      point4%x=x1
      point4%y=y1
      point4%z=z2
    
! apply the transformation to each of the points
      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point2,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point4,problem_surfaces(surface_number)%trans)
	  
! create surface triangles    
! set triangles to make the normal point outwards	  
      triangle_count=0
	  
      triangle_count=triangle_count+1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point3
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point2

      triangle_count=triangle_count+1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point4
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3

  RETURN

9000 CALL write_line('Error in build_surface_xplane:',0,.TRUE.)
     CALL write_line('x1.NE.x2',0,.TRUE.)
     CALL write_line_real('x1=',x1,0,.TRUE.)
     CALL write_line_real('x2=',x2,0,.TRUE.)
     STOP


END SUBROUTINE build_surface_xplane
