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
! SUBROUTINE read_mesh_dimensions_in_cells
!
! NAME
!     read_mesh_dimensions_in_cells
!
! DESCRIPTION
!     read mesh dimensions 
!     mesh dimension may be defined in terms of cells or real coordinates
!
! Example packet:
!
!Mesh_dimensions_in_cells
! nx ny nz                   (3*integer)
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE read_mesh_dimensions_in_cells

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: Read_mesh_dimensions_in_cells',0,output_to_screen_flag)

  read(input_file_unit,*,err=9000) nx,ny,nz

! check: nx, ny and nz should all be greater than or equal to 1
  if ( (nx.lt.1).OR.(ny.lt.1).OR.(nz.lt.1) ) GOTO 9010

  CALL write_line('FINISHED: Read_mesh_dimensions_in_cells',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error reading mesh_dimensions_in_cells packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading mesh_dimensions_in_cells packet',0,.TRUE.)
     CALL write_line('nx, ny and nz should all be greater than or equal to 1',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_mesh_dimensions_in_cells
