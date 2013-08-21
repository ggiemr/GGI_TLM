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
!     SUBROUTINE write_volume_list_vtk
!
! DESCRIPTION
!     write_volume_list_vtk:
!
!     Write tet list to vtk format file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE write_volume_list_vtk(file_unit,number_of_tets,tet_list)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: file_unit
integer	:: number_of_tets

type(xyz_tet)	:: tet_list(1:number_of_tets)

! local variables

integer 	:: tet_number
integer 	:: number_of_points
integer 	:: point_number
integer		:: vertex_number

! START

  CALL write_line('CALLED: Write_volume_list',0,output_to_screen_flag)

  number_of_points=4*number_of_tets

  write(file_unit,'(A)')'# vtk DataFile Version 2.0'
  write(file_unit,'(A)')'trim(frame_filename)'
  write(file_unit,'(A)')'ASCII'
  write(file_unit,'(A)')'DATASET POLYDATA'
  write(file_unit,'(A,I10,A)')'POINTS',number_of_points,' float'

! write point data 

    point_number=0
    
    do tet_number=1,number_of_tets
    
      do vertex_number=1,4
      
        point_number=point_number+1
	
        write(file_unit,8000)        tet_list(tet_number)%vertex(vertex_number)%x,	&
                                     tet_list(tet_number)%vertex(vertex_number)%y,	&
                                     tet_list(tet_number)%vertex(vertex_number)%z
				   
      end do ! next tet vertex   
                             
    end do! next tet 
      
8000  format(3E14.5)
  
! write tet data
    write(file_unit,'(A,2I10)')'POLYGONS',number_of_tets*4,number_of_tets*16

    point_number=0
    do tet_number=1,number_of_tets

! write 4 triangle surfaces for each tet    
      write(file_unit,8010)3,point_number  ,point_number+2,point_number+3
      write(file_unit,8010)3,point_number  ,point_number+1,point_number+2
      write(file_unit,8010)3,point_number  ,point_number+3,point_number+1
      write(file_unit,8010)3,point_number+1,point_number+3,point_number+2
      point_number=point_number+4

8010  format(I3,4I8)
      
    end do ! next tet

  CALL write_line('FINISHED: Write_volume_list',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE write_volume_list_vtk
