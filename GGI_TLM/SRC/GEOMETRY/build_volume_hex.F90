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
!SUBROUTINE build_volume_hex
!
! NAME
!     SUBROUTINE build_volume_hex
!
! DESCRIPTION
!     build_volume_hex:
!
!     create a hex volume
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/06/2025 CJS
!
!
SUBROUTINE build_volume_hex(volume_number)

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

type(xyz)	:: point1,point2,point3,point4,point5,point6,point7,point8

! START

! CREATE A TRIANGULATED RECTANGULAR BLOCK DEFINED BY OPPOSITE CORNER COORDINATES
    
      number_of_tets=5
      problem_volumes(volume_number)%number_of_tets=number_of_tets    
      ALLOCATE( problem_volumes(volume_number)%tet_list(1:number_of_tets) )

! create hex vertices
      
      point1%x=problem_volumes(volume_number)%volume_parameters(1)
      point1%y=problem_volumes(volume_number)%volume_parameters(2)
      point1%z=problem_volumes(volume_number)%volume_parameters(3)
      point2%x=problem_volumes(volume_number)%volume_parameters(4)
      point2%y=problem_volumes(volume_number)%volume_parameters(5)
      point2%z=problem_volumes(volume_number)%volume_parameters(6)
      point3%x=problem_volumes(volume_number)%volume_parameters(7)
      point3%y=problem_volumes(volume_number)%volume_parameters(8)
      point3%z=problem_volumes(volume_number)%volume_parameters(9)
      point4%x=problem_volumes(volume_number)%volume_parameters(10)
      point4%y=problem_volumes(volume_number)%volume_parameters(11)
      point4%z=problem_volumes(volume_number)%volume_parameters(12)
      
      point5%x=problem_volumes(volume_number)%volume_parameters(13)
      point5%y=problem_volumes(volume_number)%volume_parameters(14)
      point5%z=problem_volumes(volume_number)%volume_parameters(15)
      point6%x=problem_volumes(volume_number)%volume_parameters(16)
      point6%y=problem_volumes(volume_number)%volume_parameters(17)
      point6%z=problem_volumes(volume_number)%volume_parameters(18)
      point7%x=problem_volumes(volume_number)%volume_parameters(19)
      point7%y=problem_volumes(volume_number)%volume_parameters(20)
      point7%z=problem_volumes(volume_number)%volume_parameters(21)
      point8%x=problem_volumes(volume_number)%volume_parameters(22)
      point8%y=problem_volumes(volume_number)%volume_parameters(23)
      point8%z=problem_volumes(volume_number)%volume_parameters(24)
        
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
! Build the hex volume from 5 tets
	  
! set tets to make the normal point outwards	  
      tet_count=0
      
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
	  

  RETURN

END SUBROUTINE build_volume_hex
