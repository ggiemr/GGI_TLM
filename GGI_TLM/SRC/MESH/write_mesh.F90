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
!SUBROUTINE write_mesh
!
! NAME
!     SUBROUTINE write_mesh
!
! DESCRIPTION
!     write_mesh:
!
!     write volume, surface, line and point meshes to file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 12/08/2012 CJS
!
!
SUBROUTINE write_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number
integer	:: cell,number_of_cells

integer	:: surface_number
integer	:: face,number_of_faces

integer	:: line_number
integer	:: segment,number_of_cell_segments

integer	:: point_number

! START

  CALL write_line('CALLED: write_mesh',0,output_to_screen_flag)

! Open mesh file

  CALL open_file(mesh_file_unit,mesh_file_extension)

! Write general mesh parameters

  write(mesh_file_unit,*)dl,' dl'
  write(mesh_file_unit,*)nx,ny,nz,' nx ny nz'
  write(mesh_file_unit,*)mesh_xmin,mesh_xmax,' mesh_xmin,mesh_xmax'
  write(mesh_file_unit,*)mesh_ymin,mesh_ymax,' mesh_ymin,mesh_ymax'
  write(mesh_file_unit,*)mesh_zmin,mesh_zmax,' mesh_zmin,mesh_zmax'
  
  write(mesh_file_unit,*)n_volumes,' n_volumes'

  do volume_number=1,n_volumes

    number_of_cells=problem_volumes(volume_number)%number_of_cells
    
    write(mesh_file_unit,*)number_of_cells,' number_of_cells in volume=',volume_number
    
    do cell=1,number_of_cells
      write(mesh_file_unit,*)	problem_volumes(volume_number)%cell_list(cell)%cell%i,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%j,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%k,     &
				problem_volumes(volume_number)%cell_list(cell)%point
   end do ! next cell
    
  end do ! next volume number
  
  write(mesh_file_unit,*)n_surfaces,' n_surfaces'

  do surface_number=1,n_surfaces

    number_of_faces=problem_surfaces(surface_number)%number_of_faces
    
    write(mesh_file_unit,*)number_of_faces,' number_of_faces in surface=',surface_number

    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_xmin,	&
                           problem_surfaces(surface_number)%mesh_xmax,' mesh_xmin,mesh_xmax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_ymin,	&
                           problem_surfaces(surface_number)%mesh_ymax,' mesh_ymin,mesh_ymax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_zmin,	&
                           problem_surfaces(surface_number)%mesh_zmax,' mesh_zmin,mesh_zmax'

    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_xmin,	&
                           problem_surfaces(surface_number)%mesh_cell_xmax,' mesh_cell_xmin,mesh_cell_xmax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_ymin,	&
                           problem_surfaces(surface_number)%mesh_cell_ymax,' mesh_cell_ymin,mesh_cell_ymax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_zmin,	&
                           problem_surfaces(surface_number)%mesh_cell_zmax,' mesh_cell_zmin,mesh_cell_zmax'
    
    do face=1,number_of_faces
      write(mesh_file_unit,*)	problem_surfaces(surface_number)%face_list(face)%cell%i,	&
				problem_surfaces(surface_number)%face_list(face)%cell%j,	&
				problem_surfaces(surface_number)%face_list(face)%cell%k,	&
				problem_surfaces(surface_number)%face_list(face)%point
    end do ! next face
    
  end do ! next surface number
  
  write(mesh_file_unit,*)n_lines,' n_lines'

  do line_number=1,n_lines

    number_of_cell_segments=problem_lines(line_number)%number_of_cell_segments
    
    write(mesh_file_unit,*)number_of_cell_segments,' number_of_cell_segments in line=',line_number
    
    do segment=1,number_of_cell_segments
      write(mesh_file_unit,*)	problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%i,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%j,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%k,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%point,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%i,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%j,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%k,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%point
    end do ! next cell segment
    
  end do ! next line number
  
  write(mesh_file_unit,*)n_points,' n_points'

  do point_number=1,n_points

    write(mesh_file_unit,*)	problem_points(point_number)%point%x,	&
				problem_points(point_number)%point%y,	&
				problem_points(point_number)%point%z
    write(mesh_file_unit,*)	problem_points(point_number)%cell%i,	&
				problem_points(point_number)%cell%j,	&
				problem_points(point_number)%cell%k
    write(mesh_file_unit,*)	problem_points(point_number)%face%cell%i,	&
				problem_points(point_number)%face%cell%j,	&
				problem_points(point_number)%face%cell%k,	&
				problem_points(point_number)%face%point
    
  end do ! next point number
  
  CALL close_file(mesh_file_unit)

  CALL write_line('FINISHED: write_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE write_mesh
