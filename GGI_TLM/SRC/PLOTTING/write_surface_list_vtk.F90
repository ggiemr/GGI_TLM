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
!     SUBROUTINE write_surface_list_vtk
!
! DESCRIPTION
!     write_surface_list_vtk:
!
!     Write triangle list to vtk format file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE write_surface_list_vtk(file_unit,number_of_triangles,triangle_list)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: file_unit
integer	:: number_of_triangles

type(xyz_triangle)	:: triangle_list(1:number_of_triangles)

! local variables

integer 	:: triangle_number
integer 	:: number_of_points
integer 	:: point_number
integer		:: vertex_number

! START

  CALL write_line('CALLED: Write_surface_list',0,output_to_screen_flag)

  number_of_points=3*number_of_triangles

  write(file_unit,'(A)')'# vtk DataFile Version 2.0'
  write(file_unit,'(A)')'trim(frame_filename)'
  write(file_unit,'(A)')'ASCII'
  write(file_unit,'(A)')'DATASET POLYDATA'
  write(file_unit,'(A,I10,A)')'POINTS',number_of_points,' float'

! write point data 

    point_number=0
    
    do triangle_number=1,number_of_triangles
    
      do vertex_number=1,3
      
        point_number=point_number+1
	
        write(file_unit,8000)triangle_list(triangle_number)%vertex(vertex_number)%x,	&
                                     triangle_list(triangle_number)%vertex(vertex_number)%y,	&
                                     triangle_list(triangle_number)%vertex(vertex_number)%z
				   
      end do ! next triangle vertex   
                             
    end do! next triangle 
      
8000  format(3E14.5)
  
! write triangle data
    write(file_unit,'(A,2I10)')'POLYGONS',number_of_triangles,number_of_triangles*4

    point_number=0
    do triangle_number=1,number_of_triangles
    
      write(file_unit,8010)3,point_number,point_number+1,point_number+2
      point_number=point_number+3

8010  format(I3,4I8)
      
    end do ! next triangle

  CALL write_line('FINISHED: Write_surface_list',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE write_surface_list_vtk
