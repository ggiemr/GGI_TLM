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
! PROGRAM GGI_TLM_model_checks
!
! NAME
!     GGI_TLM_model_checks
!
! DESCRIPTION
!     Check meshed model before running the TLM solver:
!     Visualisation of mesh 
!     Checks intersection of different pieces of geometry 
!     Allows editing of mesh
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/08/2012 CJS
!
!
PROGRAM GGI_TLM_model_checks

USE TLM_general
USE geometry

IMPLICIT NONE

! local variables

  logical :: parallel_mesh

! START
  
  CALL write_progress('STARTED: GGI_TLM_model_checks')

  CALL write_line('TLM_model_checks',0,output_to_screen_flag)
  
  CALL write_license()
  
  CALL read_problem_name()
  
  CALL read_mesh()
   
!  CALL trim_mesh() ! redundant with new mesh generation
  
  CALL plot_mesh_volumes()
  
  CALL plot_mesh_surfaces()
  
  CALL plot_mesh_lines()
  
  CALL plot_mesh_points()
  
  CALL plot_mesh_boundary()

  CALL deallocate_geometry()
  
  CALL write_progress('FINISHED: GGI_TLM_model_checks')
  
  CALL write_line('FINISHED: GGI_TLM_model_checks',0,output_to_screen_flag)
  
END PROGRAM GGI_TLM_model_checks
