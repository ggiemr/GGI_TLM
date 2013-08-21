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
!SUBROUTINE plot_volume_material_cells
!
! NAME
!     SUBROUTINE plot_volume_material_cells
!
! DESCRIPTION
!     plot_volume_material_cells:
!
!     All the cell volumes are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/11/2012 CJS
!     revised 19/2/2013 CJS plot cells used in TLM solver only... Some may get overwritten by
!                           other materials
!
SUBROUTINE plot_volume_material_cells()

USE TLM_general
USE mesh
USE geometry_types
USE geometry
USE TLM_volume_materials
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer :: volume_material_number
integer	:: number_of_volumes,volume_number
integer	:: total_number_of_volumes
integer	:: total_number_of_cells

type(cell_point),allocatable	:: local_cell_list(:)

integer :: i,count,cell
integer :: cx,cy,cz

! START

  CALL write_line('CALLED: plot_volume_material_cells',0,output_to_screen_flag)

! read the volume material number to view  
  write(*,*)
  write(*,*)'Enter the volume material number to view'
  read(*,*)volume_material_number
  
  if ( (volume_material_number.le.0).OR.(volume_material_number.gt.n_volume_materials) ) then
    write(*,*)'volume_material_number is outside the available range'
    write(*,*)'Number of volume materials=',n_volume_materials
    RETURN
  end if
  
! work out the number of cells to write by looping over the volume list
  
  number_of_volumes=volume_material_list(volume_material_number)%n_volumes
  
  total_number_of_cells=0
  
  do cx=1,nx
    do cy=1,ny
      do cz=1,nz
	if (abs(local_cell_material(cx,cy,cz)).EQ.volume_material_number) then
	  total_number_of_cells=total_number_of_cells+1
	end if
      end do
    end do
  end do
  
  if (total_number_of_cells.EQ.0) then
    write(*,*)'Total number of cells of this material is 0'
    RETURN
  end if
  
! Allocate data to construct a local_cell_list

  ALLOCATE( local_cell_list(1:total_number_of_cells ) )
  
! fill the local_cell_list  
  
  count=0
  do cx=1,nx
    do cy=1,ny
      do cz=1,nz
	if (abs(local_cell_material(cx,cy,cz)).EQ.volume_material_number) then
	  count=count+1
	  local_cell_list(count)%cell%i=cx
	  local_cell_list(count)%cell%j=cy
	  local_cell_list(count)%cell%k=cz
	  local_cell_list(count)%point=centre
	end if
      end do
    end do
  end do
      
! open and write volume mesh to vtk format file
  CALL open_vtk_file(volume_material_cells_file_unit,volume_material_cells_file_extension,volume_material_number) 
      
  CALL write_volume_mesh_list_vtk(volume_material_cells_file_unit,	&
                                  total_number_of_cells,local_cell_list)
      
  CALL close_vtk_file(volume_material_cells_file_unit) 

  
  DEALLOCATE( local_cell_list )

  CALL write_line('FINISHED: plot_volume_material_cells',0,output_to_screen_flag)
  
  RETURN
 
  
END SUBROUTINE plot_volume_material_cells

!SUBROUTINE plot_volume_material_cells_old
!
! NAME
!     SUBROUTINE plot_volume_material_cells
!
! DESCRIPTION
!     plot_volume_material_cells:
!
!     All the cell volumes are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/11/2012 CJS
!
!
SUBROUTINE plot_volume_material_cells_old()

USE TLM_general
USE geometry_types
USE geometry
USE TLM_volume_materials
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer :: volume_material_number
integer	:: number_of_volumes,volume_number
integer	:: total_number_of_volumes
integer	:: total_number_of_cells

type(cell_point),allocatable	:: local_cell_list(:)

integer :: i,count,cell

! START

  CALL write_line('CALLED: plot_volume_material_cells',0,output_to_screen_flag)

! read the volume material number to view  
  write(*,*)
  write(*,*)'Enter the volume material number to view'
  read(*,*)volume_material_number
  
  if ( (volume_material_number.le.0).OR.(volume_material_number.gt.n_volume_materials) ) then
    write(*,*)'volume_material_number is outside the available range'
    write(*,*)'Number of volume materials=',n_volume_materials
    RETURN
  end if
  
! work out the number of cells to write by looping over the volume list
  
  number_of_volumes=volume_material_list(volume_material_number)%n_volumes
  
  total_number_of_cells=0
  
  do i=1,number_of_volumes
    volume_number=volume_material_list(volume_material_number)%volume_list(i)
    total_number_of_cells=total_number_of_cells+problem_volumes(volume_number)%number_of_cells
  end do
  
  if (total_number_of_cells.EQ.0) then
    write(*,*)'Total number of cells of this material is 0'
    RETURN
  end if
  
! Allocate data to construct a local_cell_list

  ALLOCATE( local_cell_list(1:total_number_of_cells ) )
  
! fill the local_cell_list  
  count=0
  do i=1,number_of_volumes
    volume_number=volume_material_list(volume_material_number)%volume_list(i)
    do cell=1,problem_volumes(volume_number)%number_of_cells
      count=count+1
      local_cell_list(count)=problem_volumes(volume_number)%cell_list(cell)
    end do
  end do
      
! open and write volume mesh to vtk format file
  CALL open_vtk_file(volume_material_cells_file_unit,volume_material_cells_file_extension,volume_material_number) 
      
  CALL write_volume_mesh_list_vtk(volume_material_cells_file_unit,	&
                                  total_number_of_cells,local_cell_list)
      
  CALL close_vtk_file(volume_material_cells_file_unit) 

  
  DEALLOCATE( local_cell_list )

  CALL write_line('FINISHED: plot_volume_material_cells',0,output_to_screen_flag)
  
  RETURN
 
  
END SUBROUTINE plot_volume_material_cells_old
