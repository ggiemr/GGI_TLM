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
!SUBROUTINE plot_mesh_boundary
!
! NAME
!     SUBROUTINE plot_mesh_boundary
!
! DESCRIPTION
!     plot_mesh_boundary:
!
!     write boundary to .vtk file for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!
!
SUBROUTINE plot_mesh_boundary()

USE TLM_general
USE File_information

IMPLICIT NONE

! local variables

integer 	:: number_of_faces
integer 	:: face_number
integer 	:: number_of_points
integer 	:: point_number

! START

  CALL write_line('CALLED: plot_mesh_boundary',0,output_to_screen_flag)

  CALL open_vtk_file(boundary_mesh_file_unit,boundary_mesh_file_extension,0) 

  number_of_faces=6
  number_of_points=4*number_of_faces

  write(boundary_mesh_file_unit,'(A)')'# vtk DataFile Version 2.0'
  write(boundary_mesh_file_unit,'(A)')'trim(frame_filename)'
  write(boundary_mesh_file_unit,'(A)')'ASCII'
  write(boundary_mesh_file_unit,'(A)')'DATASET POLYDATA'
  write(boundary_mesh_file_unit,'(A,I10,A)')'POINTS',number_of_points,' float'

! write point data               
! xmin face  
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymin,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymax,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymax,mesh_zmax
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymin,mesh_zmax
! xmax face  
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymin,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymax,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymax,mesh_zmax
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymin,mesh_zmax
          
! ymin face  
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymin,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymin,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymin,mesh_zmax
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymin,mesh_zmax
! ymax face  
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymax,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymax,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymax,mesh_zmax
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymax,mesh_zmax
          
! zmin face  
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymin,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymin,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymax,mesh_zmin
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymax,mesh_zmin
! zmax face  
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymin,mesh_zmax
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymin,mesh_zmax
  write(boundary_mesh_file_unit,8000)mesh_xmax,mesh_ymax,mesh_zmax
  write(boundary_mesh_file_unit,8000)mesh_xmin,mesh_ymax,mesh_zmax
      
8000 format(3E14.5)
  
! write face data
  write(boundary_mesh_file_unit,'(A,2I10)')'POLYGONS',number_of_faces,number_of_faces*5

  point_number=0
  do face_number=1,number_of_faces
    
    write(boundary_mesh_file_unit,8010)4,point_number,point_number+1,point_number+2,point_number+3
    point_number=point_number+4

8010  format(I3,4I8)
      
  end do ! next face
      
  CALL close_vtk_file(boundary_mesh_file_unit) 
    

  CALL write_line('FINISHED: plot_mesh_boundary',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_mesh_boundary
