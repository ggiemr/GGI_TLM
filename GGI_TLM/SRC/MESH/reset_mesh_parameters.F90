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
! SUBROUTINE reset_mesh_parameters
!
! NAME
!     reset_mesh_parameters
!
! DESCRIPTION
!     
!     Reset the general mesh parameters before reading the input file
!     1. mesh dimension in cells, nx,ny,nz
!     2. mesh limits mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymax,mesh_zmin,mesh_zmax
!     3. cell dimension, dl
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE reset_mesh_parameters()

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

! START  
  
  CALL write_line('CALLED: reset_mesh_parameters',0,output_to_screen_flag)

  nx=0
  ny=0
  nz=0
  
  dl=0d0
  
  mesh_xmin=0d0
  mesh_xmax=0d0
  mesh_ymin=0d0
  mesh_ymax=0d0
  mesh_zmin=0d0
  mesh_zmax=0d0
  
  CALL write_line('FINISHED: reset_mesh_parameters',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE reset_mesh_parameters
