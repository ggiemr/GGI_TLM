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
!SUBROUTINE build_surface_split_ring
!
! NAME
!     SUBROUTINE build_surface_split_ring
!
! DESCRIPTION
!     build_surface_split_ring:
!
!     create a triangulated annular split ring surface
!     
! COMMENTS
!     used to build split ring resonators
!
! HISTORY
!
!     started 16/12/2014 CJS
!
!
SUBROUTINE build_surface_split_ring(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: phi1
real*8	:: phi2
real*8	:: phi_min,phi_max,dphi
real*8	:: inner_radius,outer_radius
real*8	:: gap_angle
real*8	:: circumference
integer	:: n_circumference
integer :: n_phi
integer :: phi_loop

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3,point4

! START
    
! CREATE A TRIANGULATED ANNULAR MESH ON THE SCALE OF DL.    
      
      inner_radius=problem_surfaces(surface_number)%surface_parameters(1)
      outer_radius=problem_surfaces(surface_number)%surface_parameters(2)
      gap_angle=problem_surfaces(surface_number)%surface_parameters(3)
      
! convert gap angle to radians
      gap_angle=gap_angle*pi/180d0

! set the length of triangle edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      phi_min=gap_angle/2d0
      phi_max=2d0*pi-gap_angle/2d0
      circumference=(phi_max-phi_min)*outer_radius
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_phi=n_circumference
      
      dphi=(phi_max-phi_min)/(n_phi-1)

! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=n_phi*2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! loop over phi creating surface triangles
      
      triangle_count=0
      
      do phi_loop=1,n_phi-1
	
	phi1=phi_min+(phi_loop-1)*dphi
	phi2=phi_min+ phi_loop   *dphi
	  
! get points 1 and 2 on the outer circumference of the split_ring
	CALL rphiz_to_xyz_point(outer_radius,phi1,0d0,point1)
	CALL rphiz_to_xyz_point(outer_radius,phi2,0d0,point2)
	
! get points 3 and 4 on the inner circumference of the split_ring
	CALL rphiz_to_xyz_point(inner_radius,phi1,0d0,point3)
	CALL rphiz_to_xyz_point(inner_radius,phi2,0d0,point4)

! apply the transformation to each of the surface points
        CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
        CALL apply_transformation(point2,problem_surfaces(surface_number)%trans)
        CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
        CALL apply_transformation(point4,problem_surfaces(surface_number)%trans)

! set the two triangles for this segment of the split_ring	  
	triangle_count=triangle_count+1
	problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point2
	problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3
		
	triangle_count=triangle_count+1
	problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point3
	problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point2
	problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point4
		
      end do ! next phi

  RETURN

END SUBROUTINE build_surface_split_ring
