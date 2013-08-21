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
!SUBROUTINE write_line_mesh_list_vtk
!
! NAME
!     SUBROUTINE write_line_mesh_list_vtk
!
! DESCRIPTION
!     write_line_mesh_list_vtk:
!
!     Write line mesh list to vtk format file, each cell segment is plotted as a small tube
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 31/08/2012 CJS
!
!
SUBROUTINE write_line_mesh_list_vtk(file_unit,n_cell_segments,cell_segment_list)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer			:: file_unit
integer			:: n_cell_segments
type(cell_segment)	:: cell_segment_list(1:n_cell_segments)

! local variables

integer 	:: segment
integer 	:: point,number_of_points
integer 	:: surface,number_of_surfaces

type(xyz)	:: point1,point2

real*8		:: x1,y1,z1,x2,y2,z2

! START

  CALL write_line('CALLED: Write_line_mesh_list',0,output_to_screen_flag)

  number_of_surfaces=8  ! number of surfaces in a tube circumference

  number_of_points=number_of_surfaces*n_cell_segments*4

! write point header
  write(file_unit,'(A)')'# vtk DataFile Version 2.0'
  write(file_unit,'(A)')'trim(frame_filename)'
  write(file_unit,'(A)')'ASCII'
  write(file_unit,'(A)')'DATASET POLYDATA'
  write(file_unit,'(A,I10,A)')'POINTS',number_of_points,' float'

! loop over cell segments  
  do segment=1,n_cell_segments
  
! line end points

    CALL get_cell_point_coordinate(cell_segment_list(segment)%segment_point(1),point1)

    CALL get_cell_point_coordinate(cell_segment_list(segment)%segment_point(2),point2)
    
    CALL write_tube_points_vtk(point1%x,point1%y,point1%z,point2%x,point2%y,point2%z,	&
                               dl/10d0,number_of_surfaces,file_unit)
  
  end do ! next cell segment
  
! write point data
  write(file_unit,'(A,2I10)')'POLYGONS',number_of_surfaces*n_cell_segments,number_of_surfaces*n_cell_segments*5

  point=0
  do surface=1,number_of_surfaces*n_cell_segments
    write(file_unit,8010)4,point,point+1,point+2,point+3
    point=point+4
  end do

8010  format(I3,4I8)

  CALL write_line('FINISHED: Write_line_mesh_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error in write_line_list:',0,.TRUE.)
     CALL write_line('line length is zero',0,.TRUE.)
     STOP
  
  
END SUBROUTINE write_line_mesh_list_vtk
