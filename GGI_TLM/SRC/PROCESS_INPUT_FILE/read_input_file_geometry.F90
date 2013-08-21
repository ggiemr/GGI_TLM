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
! SUBROUTINE read_input_file_geometry
!
! NAME
!     read input file_geometry
!
! DESCRIPTION
!     read input file packets relating to the geometry and mesh
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
SUBROUTINE read_input_file_geometry

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

character*256	:: input_line

! START
    
  rewind(input_file_unit)
  
10  CONTINUE
    
! read line from input file
    read(input_file_unit,'(A)',end=100)input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,256)

    if (input_line.EQ.'mesh_outer_boundary_dimension') then
    
      CALL read_mesh_outer_boundary_dimension()
     
    else if (input_line.EQ.'mesh_dimensions_in_cells') then
    
      CALL read_mesh_dimensions_in_cells()
     
    else if (input_line.EQ.'mesh_cell_dimension') then
    
      CALL read_mesh_cell_dimension()

    else if (input_line.EQ.'volume_list') then
    
      CALL read_volume_list()

    else if (input_line.EQ.'surface_list') then
    
      CALL read_surface_list()

    else if (input_line.EQ.'line_list') then
    
      CALL read_line_list()

    else if (input_line.EQ.'point_list') then
    
      CALL read_point_list()
 
    else if (input_line.EQ.'new_mesh_generation') then
    
      new_mesh_generation=.TRUE.
      
    end if
    
    GOTO 10
    
100 CONTINUE
  
  CALL write_line('FINISHED: read_input_file_geometry',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE read_input_file_geometry
