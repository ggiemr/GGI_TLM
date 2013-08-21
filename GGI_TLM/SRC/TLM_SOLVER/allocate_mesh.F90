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
! SUBROUTINE allocate_temporary_mesh_arrays
! SUBROUTINE allocate_mesh
! SUBROUTINE deallocate_temporary_mesh_arrays
! SUBROUTINE deallocate_mesh
!
! NAME
!     allocate_temporary_mesh_arrays
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/08/2012 CJS
!
!
SUBROUTINE allocate_temporary_mesh_arrays

USE TLM_general
USE mesh
USE TLM_excitation

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: allocate_temporary_mesh_arrays',0,output_to_screen_flag)

! temporary mesh arrays to hold material codes  
  ALLOCATE ( local_surface_material(1:nx,1:ny,nzmin:nzmax,1:3) )
  local_surface_material(1:nx,1:ny,nzmin:nzmax,1:3)=0
  
  ALLOCATE ( local_cell_material(1:nx,1:ny,nzmin:nzmax) )
  local_cell_material(1:nx,1:ny,nzmin:nzmax)=0
  
  
! temporary mesh arrays to hold output codes  
  ALLOCATE ( local_surface_output(1:nx,1:ny,nzmin:nzmax,1:3) )
  local_surface_output(1:nx,1:ny,nzmin:nzmax,1:3)=0
  
  ALLOCATE ( local_cell_output(1:nx,1:ny,nzmin:nzmax) )
  local_cell_output(1:nx,1:ny,nzmin:nzmax)=0
  
  
! temporary mesh arrays to hold cable bundle codes  
  ALLOCATE ( local_surface_cable(1:nx,1:ny,nzmin:nzmax,1:3) )
  local_surface_cable(1:nx,1:ny,nzmin:nzmax,1:3)=0
 
  ALLOCATE ( local_cell_cable(1:nx,1:ny,nzmin:nzmax) )
  local_cell_cable(1:nx,1:ny,nzmin:nzmax)=0
 
 
! temporary mesh arrays to hold excitation codes  
  ALLOCATE (  local_surface_excitation(1:nx,1:ny,nzmin:nzmax,1:3) )
  local_surface_excitation(1:nx,1:ny,nzmin:nzmax,1:3)=0
  
  ALLOCATE (  local_cell_excitation(1:nx,1:ny,nzmin:nzmax) )
  local_cell_excitation(1:nx,1:ny,nzmin:nzmax)=0
  
  CALL write_line('FINISHED: allocate_temporary_mesh_arrays',0,output_to_screen_flag)

  RETURN

END SUBROUTINE allocate_temporary_mesh_arrays
!
! NAME
!     allocate_mesh
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE allocate_mesh

USE TLM_general
USE mesh
USE TLM_excitation

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: allocate_mesh',0,output_to_screen_flag)
 
  ALLOCATE ( V(1:12,0:nx+1,0:ny+1,nzmin:nzmax) )   ! TLM voltage array, link lines only
  V(1:12,1:nx,1:ny,nzmin:nzmax)=0d0

! Allocate parallel data passing variables  

  ALLOCATE(Vi_zmin(1:2*nx*ny))
  ALLOCATE(Vi_zmax(1:2*nx*ny))
  ALLOCATE(Vr_zmin(1:2*nx*ny))
  ALLOCATE(Vr_zmax(1:2*nx*ny))
  
  CALL write_line('FINISHED: allocate_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE allocate_mesh
!
! NAME
!     deallocate_temporary_mesh_arrays
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/08/2012 CJS
!
!
SUBROUTINE deallocate_temporary_mesh_arrays

USE TLM_general
USE mesh
USE TLM_excitation

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: deallocate_temporary_mesh_arrays',0,output_to_screen_flag)

! temporary mesh arrays to hold material codes  
  DEALLOCATE ( local_surface_material )
  DEALLOCATE ( local_cell_material )
  
! temporary mesh arrays to hold output codes  
  DEALLOCATE ( local_surface_output )
  DEALLOCATE ( local_cell_output )
  
! temporary mesh arrays to hold cable bundle codes  
  DEALLOCATE ( local_surface_cable )
  DEALLOCATE ( local_cell_cable )
  
! temporary mesh arrays to hold excitation codes  
  DEALLOCATE ( local_surface_excitation )
  DEALLOCATE ( local_cell_excitation )
  
  CALL write_line('FINISHED: deallocate_temporary_mesh_arrays',0,output_to_screen_flag)


  RETURN

END SUBROUTINE deallocate_temporary_mesh_arrays
