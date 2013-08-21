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
!SUBROUTINE plot_mesh_volumes
!
! NAME
!     SUBROUTINE plot_mesh_volumes
!
! DESCRIPTION
!     plot_mesh_volumes:
!
!     All the cell volumes are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
SUBROUTINE plot_mesh_volumes()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number

integer	:: number_of_cells

! START

  CALL write_line('CALLED: plot_mesh_volumes',0,output_to_screen_flag)

  do volume_number=1,n_volumes

    number_of_cells=problem_volumes(volume_number)%number_of_cells

    if (number_of_cells.gt.0) then

! open and write volume mesh to vtk format file
      CALL open_vtk_file(volume_mesh_file_unit,volume_mesh_file_extension,volume_number) 
      
      CALL write_volume_mesh_list_vtk(volume_mesh_file_unit,	&
                                  number_of_cells,problem_volumes(volume_number)%cell_list)
      
      CALL close_vtk_file(volume_mesh_file_unit) 
    
    end if ! number_of_cells.gt.0

  end do ! next volume number

  CALL write_line('FINISHED: plot_mesh_volumes',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_mesh_volumes
