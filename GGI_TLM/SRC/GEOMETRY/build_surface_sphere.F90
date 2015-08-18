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
!SUBROUTINE build_surface_sphere
!SUBROUTINE get_point_from_list
!SUBROUTINE build_surface_sphere_OLD
!
! NAME
!     SUBROUTINE build_surface_sphere
!
! DESCRIPTION
!     build_surface_sphere:
!
!     create a triangulated sphere surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!     surface roughness 13/8/2015 CJS
!
!
SUBROUTINE build_surface_sphere(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: theta,phi
real*8	:: dtheta,dphi
real*8	:: radius0,radius
real*8	:: circumference
integer	:: n_circumference
integer :: n_theta,n_phi
integer :: theta_loop,phi_loop

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3,point4

integer		:: tot_n_points
integer		:: point_count
type(xyz),allocatable	:: point_list(:)
real*8		:: r_random
integer		:: tvalue
integer		:: pvalue

! START

    
! CREATE A TRIANGULATED SPHERE MESH ON THE SCALE OF DL.    
      
      radius0=problem_surfaces(surface_number)%surface_parameters(1)

! set the length of triangle edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      circumference=2d0*pi*radius0
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_theta=n_circumference/2
      n_phi=n_circumference
      
      dtheta=pi/n_theta
      dphi=2d0*pi/n_phi
      
! calculate the number of points and allocate memory for the point data
      tot_n_points=1  				! north pole point 
      tot_n_points=tot_n_points+(n_theta-1)*n_phi  	! surface points
      tot_n_points=tot_n_points+1  		! south pole point
      
      ALLOCATE( point_list(1:tot_n_points) )

! loop over theta and phi creating the points required for the triangulation

! north pole point      
      point_count=1
      theta=0d0
      phi=0d0
      
      radius=radius0
      if (problem_surfaces(surface_number)%roughness_flag) then
        CALL random_number(r_random)
	r_random=(r_random-0.5d0)*2d0
        radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
      end if
      CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      point_list(point_count)=point1
     
      do theta_loop=2,n_theta
      
        do phi_loop=1,n_phi
	
	  point_count=point_count+1
	
	  theta=(theta_loop-1)*dtheta	  
	  phi=(phi_loop-1)*dphi
      
          radius=radius0
          if (problem_surfaces(surface_number)%roughness_flag) then
            CALL random_number(r_random)
	    r_random=(r_random-0.5d0)*2d0
            radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
          end if
          CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
          CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
          point_list(point_count)=point1
	
	end do ! next phi
	
      end do ! next theta
	
! south pole point      
      point_count=point_count+1

      theta=pi
      phi=0d0
     
      radius=radius0
      if (problem_surfaces(surface_number)%roughness_flag) then
     	CALL random_number(r_random)
	r_random=(r_random-0.5d0)*2d0
     	radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
      end if
      CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      point_list(point_count)=point1

      if (point_count.NE.tot_n_points) then
        write(*,*)'Error in build_surface_sphere'
	write(*,*)point_count,tot_n_points
	STOP
      end if
      
! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=n_theta*n_phi*2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! loop over theta and phi creating surface triangles
      
      triangle_count=0
      
      do theta_loop=1,n_theta
      
        do phi_loop=1,n_phi
	
	  CALL get_point_from_list(theta_loop  ,phi_loop  ,n_theta,n_phi,point1,point_list,tot_n_points)
	  CALL get_point_from_list(theta_loop  ,phi_loop+1,n_theta,n_phi,point2,point_list,tot_n_points)
	  CALL get_point_from_list(theta_loop+1,phi_loop+1,n_theta,n_phi,point3,point_list,tot_n_points)
	  CALL get_point_from_list(theta_loop+1,phi_loop  ,n_theta,n_phi,point4,point_list,tot_n_points)

! set triangles to make the normal point outwards	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point3
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point2
	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point4
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3
	
	end do ! next phi
	
      end do ! next theta
      
      DEALLOCATE( point_list )

  RETURN

END SUBROUTINE build_surface_sphere
!
! NAME
!     SUBROUTINE get_point_from_list
!
! DESCRIPTION
!    get_point_from_list
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/8/2015 CJS
!
!
SUBROUTINE get_point_from_list(theta,phi,n_theta,n_phi,point,point_list,tot_n_points)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer :: theta,phi,n_theta,n_phi

type(xyz)	:: point

integer		:: tot_n_points
type(xyz)	:: point_list(1:tot_n_points)


! local variables

integer :: lphi,lpoint

! START
	
  if (theta.eq.1) then
  
! point is at the north pole
    point=point_list(1)
    
  else if (theta.eq.n_theta+1) then
  
! point is at the south pole
    point=point_list(tot_n_points)
    
  else
  
    if (phi.EQ.n_phi+1) then	    
      lphi=1	      
    else 	    
      lphi=phi	      
    end if
    
    lpoint=2+(theta-2)*n_phi+(lphi-1)
    point=point_list(lpoint)
    
  end if
  
END SUBROUTINE get_point_from_list

!
!
! NAME
!     SUBROUTINE build_surface_sphere_OLD
!
! DESCRIPTION
!     build_surface_sphere:
!
!     create a triangulated sphere surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_surface_sphere_OLD(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: theta1,phi1
real*8	:: theta2,phi2
real*8	:: dtheta,dphi
real*8	:: radius
real*8	:: circumference
integer	:: n_circumference
integer :: n_theta,n_phi
integer :: theta_loop,phi_loop

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3,point4
type(xyz)	:: point5,point6,point7,point8

! START

    
! CREATE A TRIANGULATED SPHERE MESH ON THE SCALE OF DL.    
      
      radius=problem_surfaces(surface_number)%surface_parameters(1)

! set the length of triangle edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      circumference=2d0*pi*radius
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_theta=n_circumference/2
      n_phi=n_circumference
      
      dtheta=pi/n_theta
      dphi=2d0*pi/n_phi

! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=n_theta*n_phi*2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! loop over theta and phi creating surface triangles
      
      triangle_count=0
      
      do theta_loop=1,n_theta
      
        do phi_loop=1,n_phi
	
	  theta1=(theta_loop-1)*dtheta
	  theta2= theta_loop   *dtheta
	  
	  phi1=(phi_loop-1)*dphi
	  phi2= phi_loop   *dphi
	  if (phi_loop.eq.1) phi1=0d0
	  if (phi_loop.eq.n_phi) phi2=0d0

! get four points at the corners of a quad patch on the sphere surface
	  
	  CALL rthetaphi_to_xyz_point(radius,theta1,phi1,point1)
	  CALL rthetaphi_to_xyz_point(radius,theta1,phi2,point2)
	  CALL rthetaphi_to_xyz_point(radius,theta2,phi2,point3)
	  CALL rthetaphi_to_xyz_point(radius,theta2,phi1,point4)

! apply the transformation to each of the surface points
          CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
          CALL apply_transformation(point2,problem_surfaces(surface_number)%trans)
          CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
          CALL apply_transformation(point4,problem_surfaces(surface_number)%trans)

! set triangles to make the normal point outwards	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point3
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point2
	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point4
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3
	
	end do ! next phi
	
      end do ! next theta

  RETURN

END SUBROUTINE build_surface_sphere_OLD
