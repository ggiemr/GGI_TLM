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
!SUBROUTINE build_volume_sphere
!
! NAME
!     SUBROUTINE build_volume_sphere
!
! DESCRIPTION
!     build_volume_sphere:
!
!     create a tetrahedral meshed sphere volume
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_volume_sphere(volume_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: volume_number

! local variables

real*8	:: theta1,phi1
real*8	:: theta2,phi2
real*8	:: dtheta,dphi
real*8	:: radius
real*8	:: circumference
integer	:: n_circumference
integer :: n_theta,n_phi
integer :: theta_loop,phi_loop

integer	:: number_of_tets
integer	:: tet_count

type(xyz)	:: point0
type(xyz)	:: point1,point2,point3,point4
type(xyz)	:: point5,point6,point7,point8

! START

    
! CREATE A TRIANGULATED SPHERE MESH ON THE SCALE OF DL.    
      
      radius=problem_volumes(volume_number)%volume_parameters(1)

! set the length of tet edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      circumference=2d0*pi*radius
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_theta=n_circumference/2
      n_phi=n_circumference
      
      dtheta=pi/n_theta
      dphi=2d0*pi/n_phi

! calculate the number of tets and allocate memory for the tet volume data      
      number_of_tets=n_theta*n_phi*2
      problem_volumes(volume_number)%number_of_tets=number_of_tets    
      allocate( problem_volumes(volume_number)%tet_list(1:number_of_tets) )

! create central point, point0

      point0%x=0d0
      point0%y=0d0
      point0%z=0d0
      CALL apply_transformation(point0,problem_volumes(volume_number)%trans)
      
! loop over theta and phi creating volume tets
      
      tet_count=0
      
      do theta_loop=1,n_theta
      
        do phi_loop=1,n_phi
	
	  theta1=(theta_loop-1)*dtheta
	  theta2= theta_loop   *dtheta
	  
	  phi1=(phi_loop-1)*dphi
	  phi2= phi_loop   *dphi
	  if (phi_loop.eq.1) phi1=0d0
	  if (phi_loop.eq.n_phi) phi2=0d0

! get four points at the corners of a quad patch on the sphere volume
	  
	  CALL rthetaphi_to_xyz_point(radius,theta1,phi1,point1)
	  CALL rthetaphi_to_xyz_point(radius,theta1,phi2,point2)
	  CALL rthetaphi_to_xyz_point(radius,theta2,phi2,point3)
	  CALL rthetaphi_to_xyz_point(radius,theta2,phi1,point4)

! apply the transformation to each of the volume points
          CALL apply_transformation(point1,problem_volumes(volume_number)%trans)
          CALL apply_transformation(point2,problem_volumes(volume_number)%trans)
          CALL apply_transformation(point3,problem_volumes(volume_number)%trans)
          CALL apply_transformation(point4,problem_volumes(volume_number)%trans)

! set tets to make the normal point outwards	  
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
	
	end do ! next phi
	
      end do ! next theta

  RETURN

END SUBROUTINE build_volume_sphere
