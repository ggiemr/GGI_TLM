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
!SUBROUTINE build_surface_annulus
!
! NAME
!     SUBROUTINE build_surface_annulus
!
! DESCRIPTION
!     build_surface_annulus:
!
!     create a triangulated annular surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/12/2014 CJS
!
!
SUBROUTINE build_surface_annulus(surface_number)

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
real*8	:: dphi
real*8	:: inner_radius,outer_radius
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

! set the length of triangle edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      circumference=2d0*pi*outer_radius
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_phi=n_circumference
      
      dphi=2d0*pi/n_phi

! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=n_phi*2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! loop over phi creating surface triangles
      
      triangle_count=0
      
      do phi_loop=1,n_phi
	
	phi1=(phi_loop-1)*dphi
	phi2= phi_loop   *dphi
	if (phi_loop.eq.1) phi1=0d0
	if (phi_loop.eq.n_phi) phi2=0d0
	  
! get points 1 and 2 on the outer circumference of the annulus
	CALL rphiz_to_xyz_point(outer_radius,phi1,0d0,point1)
	CALL rphiz_to_xyz_point(outer_radius,phi2,0d0,point2)
	
! get points 3 and 4 on the inner circumference of the annulus
	CALL rphiz_to_xyz_point(inner_radius,phi1,0d0,point3)
	CALL rphiz_to_xyz_point(inner_radius,phi2,0d0,point4)

! apply the transformation to each of the surface points
        CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
        CALL apply_transformation(point2,problem_surfaces(surface_number)%trans)
        CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
        CALL apply_transformation(point4,problem_surfaces(surface_number)%trans)

! set the two triangles for this segment of the annulus	  
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

END SUBROUTINE build_surface_annulus
