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
!SUBROUTINE mesh_tet
!
! NAME
!     SUBROUTINE mesh_tet
!
! DESCRIPTION
!     mesh_tet:
!
!     The new algorithm is based on meshing the vertices, then the four triangular faces,
!     then volume filling inside the surfaces to give the cells
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!     rewritten 14/5/2013 CJS - new strategy which provides a much more robust 
!                               volume mesh, it is also consistent with the triangle mesh generation strategy
!                               so that volumes with surface coatings are meshed consistently.
!
!
  SUBROUTINE mesh_tet_new(volume_mesh,tet)
  
USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE file_information
USE local_mesh
USE constants

IMPLICIT NONE

  integer		:: volume_mesh(1:nx,1:ny,1:nz)
  type(xyz_tet)		:: tet

! local variables

  integer 		:: ix,iy,iz
  
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
  real*8		:: x4,y4,z4
 
  type(xyz)		:: min_point,max_point
  
  type(cell_point)	:: mesh_cell
    
  logical		:: inside_mesh
    
! START

!  CALL write_line('CALLED: mesh_tet',0,output_to_screen_flag)
  
! 1. Extract the vertices  
  
  x1=tet%vertex(1)%x
  y1=tet%vertex(1)%y
  z1=tet%vertex(1)%z
  
  x2=tet%vertex(2)%x
  y2=tet%vertex(2)%y
  z2=tet%vertex(2)%z
  
  x3=tet%vertex(3)%x
  y3=tet%vertex(3)%y
  z3=tet%vertex(3)%z
  
  x4=tet%vertex(4)%x
  y4=tet%vertex(4)%y
  z4=tet%vertex(4)%z
  
! 2. Get the extent of the tet in x,y and z

  min_point%x=min(x1,x2,x3,x4)
  min_point%y=min(y1,y2,y3,y4)
  min_point%z=min(z1,z2,z3,z4)
  
  max_point%x=max(x1,x2,x3,x4)
  max_point%y=max(y1,y2,y3,y4)
  max_point%z=max(z1,z2,z3,z4)
     
! 4. Allocate a local grid for the triangle then reset the grid data

  CALL allocate_local_grid(min_point,max_point)

  CALL allocate_local_grid_tet(min_point,max_point)
  
! 5. mesh the tet cells

  CALL local_grid_tetrahedron(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)

! 6. Transfer the tet cells into the volume_mesh array

! go throught the local grid transferring the cells

  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
      
        if (local_grid_tet(ix,iy,iz,local_centre).NE.0) then  
	
  	  mesh_cell%cell%i=ix
  	  mesh_cell%cell%j=iy
  	  mesh_cell%cell%k=iz
  	  mesh_cell%point=centre
	
          CALL cell_point_inside_mesh(mesh_cell,inside_mesh)
	
	  if (inside_mesh) then ! set face in the volume mesh array
	    volume_mesh(mesh_cell%cell%i,mesh_cell%cell%j,mesh_cell%cell%k)=1
	  end if
	  
	end if ! this cell is set

      end do
    end do
  end do

! 7. Deallocate local grid

  DEALLOCATE( local_grid )

  DEALLOCATE( local_grid_tet )
  
!  CALL write_line('FINISHED: mesh_tet',0,output_to_screen_flag)

  RETURN
  
  END SUBROUTINE mesh_tet_new
