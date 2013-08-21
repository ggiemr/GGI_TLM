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
!SUBROUTINE build_volume_pyramid_ram
!
! NAME
!     SUBROUTINE build_volume_pyramid_ram
!
! DESCRIPTION
!     build_volume_pyramid_ram:
!
!     create a pyramid_ram volume from tets
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/05/2013 CJS
!
!
SUBROUTINE build_volume_pyramid_ram(volume_number)

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

real*8  :: height,half_side_length,base_height

type(xyz)	:: point1,point2,point3,point4,point5,point6,point7,point8,point9

! START

! CREATE A TRIANGULATED RECTANGULAR BLOCK DEFINED BY OPPOSITE CORNER COORDINATES
    
      number_of_tets=7
      problem_volumes(volume_number)%number_of_tets=number_of_tets    
      ALLOCATE( problem_volumes(volume_number)%tet_list(1:number_of_tets) )

! create pyramid vertices

      half_side_length=problem_volumes(volume_number)%volume_parameters(1)
      height=problem_volumes(volume_number)%volume_parameters(2)
      base_height=problem_volumes(volume_number)%volume_parameters(3)

! square base
      point1%x=-half_side_length
      point1%y=-half_side_length
      point1%z=0d0
      
      point2%x=+half_side_length
      point2%y=-half_side_length
      point2%z=0d0
      
      point3%x=+half_side_length
      point3%y=+half_side_length
      point3%z=0d0
      
      point4%x=-half_side_length
      point4%y=+half_side_length
      point4%z=0d0
      
! square top of base
      point5%x=-half_side_length
      point5%y=-half_side_length
      point5%z=base_height
      
      point6%x=+half_side_length
      point6%y=-half_side_length
      point6%z=base_height
      
      point7%x=+half_side_length
      point7%y=+half_side_length
      point7%z=base_height
      
      point8%x=-half_side_length
      point8%y=+half_side_length
      point8%z=base_height
      
      point9%x=0d0
      point9%y=0d0
      point9%z=height
    
! apply the transformation to each of the points
      CALL apply_transformation(point1,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point2,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point3,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point4,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point5,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point6,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point7,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point8,problem_volumes(volume_number)%trans)
      CALL apply_transformation(point9,problem_volumes(volume_number)%trans)
	  
! create volume tets    
! set tets to make the normal point outwards	  
      tet_count=0

! five tets defining the square base 
      
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point5    
      
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point3
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point7    
      
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point7    
      
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point6
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point7    
      
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point4
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point8
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point7    
	  
! two tets creating the pyramid point      
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point6
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point9    
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point5
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point7
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point8
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point9    

  RETURN

END SUBROUTINE build_volume_pyramid_ram
