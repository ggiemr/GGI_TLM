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
!SUBROUTINE plot_mesh_points
!
! NAME
!     SUBROUTINE plot_mesh_points
!
! DESCRIPTION
!     plot_mesh_points:
!
!     All the triangulated points are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE plot_mesh_points()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: point_number

! START

  CALL write_line('CALLED: plot_mesh_points',0,output_to_screen_flag)

  do point_number=1,n_points

! open and write point to vtk format file
    CALL open_vtk_file(point_mesh_file_unit,point_mesh_file_extension,point_number) 
    
    CALL write_point_mesh_list_vtk(point_mesh_file_unit,problem_points(point_number)%cell)
    
    CALL close_vtk_file(point_mesh_file_unit) 

  end do ! next point number

  CALL write_line('FINISHED: plot_mesh_points',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_mesh_points
