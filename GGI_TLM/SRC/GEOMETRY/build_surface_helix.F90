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
!SUBROUTINE build_surface_helix
!
! NAME
!     SUBROUTINE build_surface_helix
!
! DESCRIPTION
!     build_surface_helix:
!
!     create a triangulated helical surface
!     
! COMMENTS
!     used to build helical antennas
!
! HISTORY
!
!     started 1/5/2018 CJS
!
!
SUBROUTINE build_surface_helix(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: height
real*8	:: phi0
real*8	:: phi1
real*8	:: phi2
real*8	:: pitch
real*8	:: dphi
real*8	:: radius
real*8	:: l_circumference
integer	:: n_circumference
integer	:: n_z
integer :: n_phi
integer :: phi_loop
integer :: z_loop
integer :: dirn
real*8  :: dz
real*8  :: dphiz
real*8	:: phia,phib,phic,phid
real*8  :: z1,z2

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3,point4

! START
    
! CREATE A TRIANGULATED HELICAL MESH ON THE SCALE OF DL.    
      
      radius=problem_surfaces(surface_number)%surface_parameters(1)
      height=problem_surfaces(surface_number)%surface_parameters(2)
      phi1=problem_surfaces(surface_number)%surface_parameters(3)
      phi2=problem_surfaces(surface_number)%surface_parameters(4)
      pitch=problem_surfaces(surface_number)%surface_parameters(5)
      
      dirn=+1
      if (problem_surfaces(surface_number)%surface_parameters(6).LT.0d0) then
        dirn=-1
      end if
      
! convert angles to radians
      phi1=phi1*pi/180d0
      phi2=phi2*pi/180d0

! set the length of triangle edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference
      
      l_circumference=(phi2-phi1)*radius
      
! calculate an approximate number of edge lengths around the circumference
      n_circumference=4*NINT(l_circumference/dl)      
      if (n_circumference.lt.4) n_circumference=4
      
! calculate an approximate number of edge lengths in z
      n_z=4*NINT(height/dl)    
      if (n_z.lt.4) n_z=4

! calculate the number of edges along the circumferential direction      
      n_phi=n_circumference

! phi step in circumferential direction   
      dphi=(phi2-phi1)/n_circumference
      
! z step      
      dz=height/n_z
      
! phi step when moving to the next z    
      dphiz=dirn*(2d0*pi*height/pitch)/(n_z-1)

! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=2*n_phi*n_z
      
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )
      
      triangle_count=0

! loop over z 
      
      do z_loop=1,n_z
        
        z1=(z_loop-1)*dz
        z2= z_loop   *dz

! calculate the initial value of phi for this z        
        phi0=phi1+(z_loop-1)*dphiz
	
! loop over phi (circumferential direction) creating surface triangles
            
        do phi_loop=1,n_phi

  	  phia=phi0+(phi_loop-1)*dphi
	  phib=phi0+(phi_loop  )*dphi
	  phic=phi0+(phi_loop-1)*dphi+dphiz
	  phid=phi0+(phi_loop  )*dphi+dphiz
	  
! get points a and b on the helix
	  CALL rphiz_to_xyz_point(radius,phia,z1,point1)
	  CALL rphiz_to_xyz_point(radius,phib,z1,point2)
	
! get points c and d on the helix
	  CALL rphiz_to_xyz_point(radius,phic,z2,point3)
	  CALL rphiz_to_xyz_point(radius,phid,z2,point4)

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
		
      end do ! next z

  RETURN

END SUBROUTINE build_surface_helix
