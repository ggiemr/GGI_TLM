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
!     SUBROUTINE write_point_list_vtk
!
! DESCRIPTION
!     write_point_list_vtk:
!
!     Write point list to vtk format file, each point is plotted as a small diamond shape
!     made up of 8 triangles
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE write_point_list_vtk(file_unit,point)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: file_unit

type(xyz)	:: point

! local variables

integer 	:: number_of_points
integer 	:: number_of_triangles

real*8		:: offset

! START

  CALL write_line('CALLED: Write_point_list',0,output_to_screen_flag)

  number_of_triangles=8

  number_of_points=6

  write(file_unit,'(A)')'# vtk DataFile Version 2.0'
  write(file_unit,'(A)')'trim(frame_filename)'
  write(file_unit,'(A)')'ASCII'
  write(file_unit,'(A)')'DATASET POLYDATA'
  write(file_unit,'(A,I10,A)')'POINTS',number_of_points,' float'

! write point data 
  
  offset=dl/5d0
	
  write(file_unit,8000)point%x		,point%y	,point%z+offset
  write(file_unit,8000)point%x+offset	,point%y	,point%z
  write(file_unit,8000)point%x		,point%y+offset	,point%z
  write(file_unit,8000)point%x-offset	,point%y	,point%z
  write(file_unit,8000)point%x		,point%y-offset	,point%z
  write(file_unit,8000)point%x		,point%y	,point%z-offset
      
8000  format(3E14.5)
  
! write triangle data
  write(file_unit,'(A,2I10)')'POLYGONS',number_of_triangles,number_of_triangles*4

  write(file_unit,8010)3,0,1,2
  write(file_unit,8010)3,0,2,3
  write(file_unit,8010)3,0,3,4
  write(file_unit,8010)3,0,4,1
  write(file_unit,8010)3,5,2,1
  write(file_unit,8010)3,5,3,2
  write(file_unit,8010)3,5,4,3
  write(file_unit,8010)3,5,1,4

8010  format(I3,4I8)

  CALL write_line('FINISHED: Write_point_list',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE write_point_list_vtk
