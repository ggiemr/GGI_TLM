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
!SUBROUTINE build_volume_cylinder
!
! NAME
!     SUBROUTINE build_volume_cylinder
!
! DESCRIPTION
!     build_volume_cylinder:
!
!     create a triangulated cylinder volume
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_volume_cylinder(volume_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: volume_number

! local variables

real*8	:: phi1
real*8	:: phi2
real*8	:: zmin,zmax
real*8	:: dphi
real*8	:: radius
real*8	:: circumference
integer	:: n_circumference
integer :: n_phi
integer :: phi_loop

integer	:: number_of_tets
integer	:: tet_count

type(xyz)	:: point0
type(xyz)	:: point1,point2,point3,point4,point5,point6

! START
    
! CREATE A TRIANGULATED CYLINDER MESH ON THE SCALE OF DL.    
      
      radius=problem_volumes(volume_number)%volume_parameters(1)
      zmin=-problem_volumes(volume_number)%volume_parameters(2)/2d0
      zmax= problem_volumes(volume_number)%volume_parameters(2)/2d0

! set the length of tet edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      circumference=2d0*pi*radius
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_phi=n_circumference
      
      dphi=2d0*pi/n_phi

! calculate the number of tets and allocate memory for the triangulated volume data      
      number_of_tets=n_phi*4
      problem_volumes(volume_number)%number_of_tets=number_of_tets    
      allocate( problem_volumes(volume_number)%tet_list(1:number_of_tets) )

! create central point, point0

      point0%x=0d0
      point0%y=0d0
      point0%z=0d0
      CALL apply_transformation(point0,problem_volumes(volume_number)%trans)

! loop over phi creating volume tets
      
      tet_count=0
      
      do phi_loop=1,n_phi
	
	phi1=(phi_loop-1)*dphi
	phi2= phi_loop   *dphi
	if (phi_loop.eq.1) phi1=0d0
	if (phi_loop.eq.n_phi) phi2=0d0

! get point1 at the centre of the top cylinder volume
	point1%x=0d0
	point1%y=0d0
	point1%z=zmax
	  
! get points 2 and 3 on the circumference of the cylinder top
	CALL rphiz_to_xyz_point(radius,phi1,zmax,point2)
	CALL rphiz_to_xyz_point(radius,phi2,zmax,point3)
	  
! get points 4 and 5 on the circumference of the cylinder bottom
	CALL rphiz_to_xyz_point(radius,phi1,zmin,point4)
	CALL rphiz_to_xyz_point(radius,phi2,zmin,point5)

! get point6 at the centre of the top cylinder volume
	point6%x=0d0
	point6%y=0d0
	point6%z=zmin

! apply the transformation to each of the volume points
        CALL apply_transformation(point1,problem_volumes(volume_number)%trans)
        CALL apply_transformation(point2,problem_volumes(volume_number)%trans)
        CALL apply_transformation(point3,problem_volumes(volume_number)%trans)
        CALL apply_transformation(point4,problem_volumes(volume_number)%trans)
        CALL apply_transformation(point5,problem_volumes(volume_number)%trans)
        CALL apply_transformation(point6,problem_volumes(volume_number)%trans)

! set tets to make the normal point outwards	  
	tet_count=tet_count+1
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point2
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point3
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	
	tet_count=tet_count+1
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point4
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point5
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	
	tet_count=tet_count+1
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point5
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point3
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	
	tet_count=tet_count+1
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point6
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point5
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point4
	problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	
      end do ! next phi

  RETURN

END SUBROUTINE build_volume_cylinder
