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
! PROGRAM GGI_TLM_model_builder
!
! NAME
!     GGI_TLM_model_builder
!
! DESCRIPTION
!     Builds the geometry based on tetrahedral volumes, triangulated surfaces, line segments and points
!     These unstructured geometric objects are then meshed using a uniform cubic mesh.
!
!     Meshed volumes are cell based, Meshed surfaces are face based, 
!     line segments are based on cell centre to cell face centre line segments (this
!     reflects the TLM cable model algorithm) and points are cell centre based.
!
!     The model builder outputs the unstructured geometry representation to vtk format files
!     suitable for visualisation with ParaView
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
PROGRAM GGI_TLM_model_builder

USE TLM_general
USE geometry

IMPLICIT NONE

! local variables

! START

  write(*,*)'GGI_TLM_model_builder'
  
  CALL write_progress('STARTED: GGI_TLM_model_builder')
  
  CALL write_line('GGI_TLM_model_builder',0,output_to_screen_flag)
  
  CALL write_license()
  
  CALL reset_packet_data()
  
  CALL reset_mesh_parameters()
  
  CALL read_problem_name()

  CALL open_general_files()
  
  CALL read_input_file_geometry()
  
  CALL set_mesh_parameters()

! BUILD GEOMETRY
  
  CALL build_volume_geometry()
  
  CALL build_surface_geometry()
  
  CALL build_line_geometry()
  
  CALL build_point_geometry()

! BUILD MESH
  
  CALL build_volume_mesh()
  
  CALL build_surface_mesh()
  
  CALL build_line_mesh()
  
  CALL build_point_mesh()

! WRITE GEOMETRY AND MESH FILES

  CALL write_mesh()

  CALL deallocate_geometry() ! deallocates geometry and mesh data

  CALL close_general_files()
  
  CALL write_progress('FINISHED: GGI_TLM_model_builder')

  CALL write_line('FINISHED: GGI_TLM_model_builder',0,output_to_screen_flag)
  
END PROGRAM GGI_TLM_model_builder
