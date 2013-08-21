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
!SUBROUTINE mesh_triangle
!
! NAME
!     SUBROUTINE mesh_triangle
!
! DESCRIPTION
!     mesh_triangle:
!
!     mesh a single triangle.
!     The new algorithm is based on meshing the vertices, then edges then filling the surfaces
!     inside the edges. 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/08/2012 CJS
!     revised 14/9/2012 CJS   - simpler strategy which avoids duplicate faces
!     rewritten 14/5/2013 CJS - new strategy which provides a much more robust 
!                               surface mesh, it is also consistent with the tet mesh generation strategy
!                               so that volumes with surface coatings are meshed consistently.
!
SUBROUTINE mesh_triangle_new(surface_mesh,triangle)
  
USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE file_information
USE local_mesh
USE constants

IMPLICIT NONE

  integer		:: surface_mesh(1:nx,1:ny,1:nz,1:6)
  type(xyz_triangle)	:: triangle

! local variables

  integer 		:: ix,iy,iz
   
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
  
  integer 		:: ix1,iy1,iz1
  integer 		:: ix2,iy2,iz2
  integer 		:: ix3,iy3,iz3
   
  type(xyz)		:: min_point,max_point
  
  type(xyz)		:: normal
  
  type(cell_point)	:: mesh_face
  
  logical		:: inside
    
! START
  
! 1. Extract the vertices  
  
  x1=triangle%vertex(1)%x
  y1=triangle%vertex(1)%y
  z1=triangle%vertex(1)%z
  
  x2=triangle%vertex(2)%x
  y2=triangle%vertex(2)%y
  z2=triangle%vertex(2)%z
  
  x3=triangle%vertex(3)%x
  y3=triangle%vertex(3)%y
  z3=triangle%vertex(3)%z
  
! 2. Get the extent of the triangle in x,y and z

  min_point%x=min(x1,x2,x3)
  min_point%y=min(y1,y2,y3)
  min_point%z=min(z1,z2,z3)
  
  max_point%x=max(x1,x2,x3)
  max_point%y=max(y1,y2,y3)
  max_point%z=max(z1,z2,z3)
     
! 4. Allocate a local grid for the triangle then reset the grid data

  CALL allocate_local_grid(min_point,max_point)
  
! 5. mesh the triangle surfaces

  CALL local_grid_triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3)

! 6. Transfer the triangle surfaces into the surface_mesh array

! get the triangle normal direction

  CALL triangle_normal(triangle,normal)

! go throught the local grid transferring the surfaces

  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
      
! surfaces normal to x
      
        if (local_grid(ix,iy,iz,local_xface).NE.0) then  
! the xmin face of this cell is set

          mesh_face%cell%i=ix
          mesh_face%cell%j=iy
          mesh_face%cell%k=iz

	  if ( (normal%x.gt.0d0) ) then
            mesh_face%cell%i=mesh_face%cell%i-1
	    mesh_face%point=face_xmax
	  else
            mesh_face%cell%i=mesh_face%cell%i
	    mesh_face%point=face_xmin
	  end if
	
	  CALL cell_point_inside_mesh(mesh_face,inside)
	
	  if (inside) then ! set face in the surface mesh array
	    surface_mesh(mesh_face%cell%i,mesh_face%cell%j,mesh_face%cell%k,mesh_face%point)=1
	  end if
	
	end if ! surface normal to x
      
! surfaces normal to y
      
        if (local_grid(ix,iy,iz,local_yface).NE.0) then  
! the ymin face of this cell is set

          mesh_face%cell%i=ix
          mesh_face%cell%j=iy
          mesh_face%cell%k=iz

	  if ( (normal%y.gt.0d0) ) then
            mesh_face%cell%j=mesh_face%cell%j-1
	    mesh_face%point=face_ymax
	  else
            mesh_face%cell%j=mesh_face%cell%j
	    mesh_face%point=face_ymin
	  end if
	
	  CALL cell_point_inside_mesh(mesh_face,inside)
	
	  if (inside) then ! set face in the surface mesh array
	    surface_mesh(mesh_face%cell%i,mesh_face%cell%j,mesh_face%cell%k,mesh_face%point)=1
	  end if
	
	end if ! surface normal to y
      
! surfaces normal to z
      
        if (local_grid(ix,iy,iz,local_zface).NE.0) then  
! the zmin face of this cell is set

          mesh_face%cell%i=ix
          mesh_face%cell%j=iy
          mesh_face%cell%k=iz

	  if ( (normal%z.gt.0d0) ) then
            mesh_face%cell%k=mesh_face%cell%k-1
	    mesh_face%point=face_zmax
	  else
            mesh_face%cell%k=mesh_face%cell%k
	    mesh_face%point=face_zmin
	  end if
	
	  CALL cell_point_inside_mesh(mesh_face,inside)
	
	  if (inside) then ! set face in the surface mesh array
	    surface_mesh(mesh_face%cell%i,mesh_face%cell%j,mesh_face%cell%k,mesh_face%point)=1
	  end if
	
	end if ! surface normal to z

      end do
    end do
  end do

! 7. Deallocate local grid

  DEALLOCATE( local_grid )
  
  RETURN
  
  END SUBROUTINE mesh_triangle_new
