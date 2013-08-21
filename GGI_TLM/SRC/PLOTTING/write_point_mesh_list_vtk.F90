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
!
! NAME
!     SUBROUTINE write_point_mesh_list_vtk
!
! DESCRIPTION
!     write_point_mesh_list_vtk:
!
!     Write face list to vtk format file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/08/2012 CJS
!
!
SUBROUTINE write_point_mesh_list_vtk(file_unit,mesh_cell)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: file_unit

type(cell_point)	:: mesh_cell

! local variables

integer 	:: number_of_points
type(xyz)	:: point

! START

  CALL write_line('CALLED: Write_point_mesh_list',0,output_to_screen_flag)

  number_of_points=8

  write(file_unit,'(A)')'# vtk DataFile Version 2.0'
  write(file_unit,'(A)')'trim(frame_filename)'
  write(file_unit,'(A)')'ASCII'
  write(file_unit,'(A)')'DATASET POLYDATA'
  write(file_unit,'(A,I10,A)')'POINTS',number_of_points,' float'

! write point data 

  CALL get_cell_centre_coordinate(mesh_cell%cell,point)    
	
  write(file_unit,8000)point%x-dl/2d0,point%y-dl/2d0,point%z-dl/2d0
  write(file_unit,8000)point%x+dl/2d0,point%y-dl/2d0,point%z-dl/2d0
  write(file_unit,8000)point%x+dl/2d0,point%y+dl/2d0,point%z-dl/2d0
  write(file_unit,8000)point%x-dl/2d0,point%y+dl/2d0,point%z-dl/2d0
  write(file_unit,8000)point%x-dl/2d0,point%y-dl/2d0,point%z+dl/2d0
  write(file_unit,8000)point%x+dl/2d0,point%y-dl/2d0,point%z+dl/2d0
  write(file_unit,8000)point%x+dl/2d0,point%y+dl/2d0,point%z+dl/2d0
  write(file_unit,8000)point%x-dl/2d0,point%y+dl/2d0,point%z+dl/2d0
      
8000  format(3E14.5)
  
! write face data
  write(file_unit,'(A,2I10)')'POLYGONS',6,30
    
  write(file_unit,8010)4,0,3,2,1
  write(file_unit,8010)4,4,5,6,7
  write(file_unit,8010)4,0,1,5,4
  write(file_unit,8010)4,2,3,7,6
  write(file_unit,8010)4,1,2,6,5
  write(file_unit,8010)4,0,4,7,3
  
8010  format(I3,4I8)

  CALL write_line('FINISHED: Write_point_mesh_list',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE write_point_mesh_list_vtk
