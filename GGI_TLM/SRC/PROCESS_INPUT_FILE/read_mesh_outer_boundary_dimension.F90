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
! SUBROUTINE read_mesh_outer_boundary_dimension
!
! NAME
!     read_mesh_outer_boundary_dimension
!
! DESCRIPTION
!     read mesh outer boundary dimension
!
! Example packet format
!
!Mesh_outer_boundary_dimension
! mesh_xmin mesh_xmax   mesh_ymin mesh_ymax   mesh_zmin mesh_zmax
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE read_mesh_outer_boundary_dimension

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: Read_mesh_outer_boundary_dimension',0,output_to_screen_flag)

  read(input_file_unit,*,err=9000)mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymax,mesh_zmin,mesh_zmax

! checks

  CALL write_line('FINISHED: Read_mesh_outer_boundary_dimension',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error reading read_mesh_outer_boundary_dimension packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_mesh_outer_boundary_dimension
