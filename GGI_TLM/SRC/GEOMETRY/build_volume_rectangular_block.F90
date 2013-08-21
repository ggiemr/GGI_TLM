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
!SUBROUTINE build_volume_rectangular_block
!
! NAME
!     SUBROUTINE build_volume_rectangular_block
!
! DESCRIPTION
!     build_volume_rectangular_block:
!
!     create a triangulated rectangular block volume
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_volume_rectangular_block(volume_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: volume_number

! local variables

integer	:: number_of_tets
integer	:: tet_count

real*8	:: lx,ly,lz

type(xyz)	:: point0
type(xyz)	:: point1,point2,point3,point4
type(xyz)	:: point5,point6,point7,point8

! START

! CREATE A TRIANGULATED RECTANGULAR BLOCK  
    
      number_of_tets=12
      problem_volumes(volume_number)%number_of_tets=number_of_tets    
      allocate( problem_volumes(volume_number)%tet_list(1:number_of_tets) )

! create rectangular_block vertices
      lx=problem_volumes(volume_number)%volume_parameters(1)
      ly=problem_volumes(volume_number)%volume_parameters(2)
      lz=problem_volumes(volume_number)%volume_parameters(3)

! create central point, point0

      point0%x=0d0
      point0%y=0d0
      point0%z=0d0
      CALL apply_transformation(point0,problem_volumes(volume_number)%trans)

      point1%x=-lx/2d0
      point1%y=-ly/2d0
      point1%z=-lz/2d0
      point2%x=+lx/2d0
      point2%y=-ly/2d0
      point2%z=-lz/2d0
      point3%x=+lx/2d0
      point3%y=+ly/2d0
      point3%z=-lz/2d0
      point4%x=-lx/2d0
      point4%y=+ly/2d0
      point4%z=-lz/2d0
      point5%x=-lx/2d0
      point5%y=-ly/2d0
      point5%z=+lz/2d0
      point6%x=+lx/2d0
      point6%y=-ly/2d0
      point6%z=+lz/2d0
      point7%x=+lx/2d0
      point7%y=+ly/2d0
      point7%z=+lz/2d0
      point8%x=-lx/2d0
      point8%y=+ly/2d0
      point8%z=+lz/2d0
    
! apply the transformation to each of the points
      CALL apply_transformation(point1,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point2,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point3,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point4,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point5,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point6,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point7,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point8,problem_volumes(volume_number)%trans)
	  
! create volume tets    
! set tets to make the normal point outwards	  
      tet_count=0
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point3
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point3
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point6
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point6
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point3
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point6
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point6
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point8
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point3
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point8
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point8
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point8
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0   

  RETURN

END SUBROUTINE build_volume_rectangular_block
