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
! MODULE File_general
! MODULE File_information
! MODULE file_header
! MODULE output_formats
!
! NAME
!     MODULE File_general
!
! DESCRIPTION
!     general file reading information, filename length and line lengths
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
MODULE File_general

IMPLICIT NONE

  integer,parameter :: filename_length=256
  
  integer,parameter :: line_length=256

END MODULE File_general

!
! NAME
!     MODULE File_information
!
! DESCRIPTION
!     file information: unit numbers and file extensions
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
MODULE File_information

IMPLICIT NONE
  
  character(len=4)	:: input_file_extension='.inp'
  integer,parameter	:: input_file_unit=10
  
  character(len=8)	:: warning_file_extension='.warning'
  integer,parameter	:: warning_file_unit=11
  
  character(len=11)	:: tet_volume_file_extension='.v_geom.vtk'
  integer,parameter	:: tet_volume_file_unit=12
  
  character(len=11)	:: volume_mesh_file_extension='.v_mesh.vtk'
  integer,parameter	:: volume_mesh_file_unit=13
  
  character(len=11)	:: triangulated_surface_file_extension='.s_geom.vtk'
  integer,parameter	:: triangulated_surface_file_unit=14
  
  character(len=11)	:: surface_mesh_file_extension='.s_mesh.vtk'
  integer,parameter	:: surface_mesh_file_unit=15
  
  character(len=11)	:: line_segment_file_extension='.l_geom.vtk'
  integer,parameter	:: line_segment_file_unit=16
  
  character(len=11)	:: line_mesh_file_extension='.l_mesh.vtk'
  integer,parameter	:: line_mesh_file_unit=17
  
  character(len=11)	:: point_file_extension='.p_geom.vtk'
  integer,parameter	:: point_file_unit=18
  
  character(len=11)	:: point_mesh_file_extension='.p_mesh.vtk'
  integer,parameter	:: point_mesh_file_unit=19
  
  character(len=11)	:: boundary_mesh_file_extension='.b_mesh.vtk'
  integer,parameter	:: boundary_mesh_file_unit=20
  
  character(len=15)	:: volume_material_cells_file_extension='.vmat_cells.vtk'
  integer,parameter	:: volume_material_cells_file_unit=21
    
  character(len=15)	:: surface_material_faces_file_extension='.smat_faces.vtk'
  integer,parameter	:: surface_material_faces_file_unit=22
  
  character(len=11)	:: cable_mesh_file_extension='.c_mesh.vtk'
  integer,parameter	:: cable_mesh_file_unit=23
  
  integer,parameter	:: mode_file_unit=24
  
  character(len=5)	:: mesh_file_extension='.mesh'
  integer,parameter	:: mesh_file_unit=30
  
  character(len=5)	:: surface_material_file_extn='.smat'
  integer,parameter	:: surface_material_file_unit=31
  
  character(len=5)	:: volume_material_file_extn='.vmat'
  integer,parameter	:: volume_material_file_unit=32
  
  character(len=6)	:: cable_geometry_file_extn='.cable'
  integer,parameter	:: cable_geometry_file_unit=33
  
  character(len=5)	:: info_file_extn='.info'
  integer,parameter	:: info_file_unit=40
  
  character(len=11)	:: cable_info_file_extn='.cable_info'
  integer,parameter	:: cable_info_file_unit=41
  
  character(len=12)	:: cable_model_file_extn='.cable_model'
  integer,parameter	:: cable_model_file_unit=42
 
  character(len=11)	:: field_output_extn='.field.tout'
  integer,parameter 	:: field_output_unit=50
  
  character(len=16)	:: excitation_output_extn='.excitation.tout'
  integer,parameter 	:: excitation_output_unit=51
 
  character(len=19)	:: cable_current_output_extn='.cable_current.tout'
  integer,parameter 	:: cable_current_output_unit=52
 
  character(len=18)	:: volume_field_output_extn='.volume_field.tout'
  integer,parameter 	:: volume_field_output_unit=53
 
  character(len=26)	:: volume_average_field_output_extn='.volume_average_field.tout'
  integer,parameter 	:: volume_average_field_output_unit=54
  
  character(len=19)	:: surface_field_output_extn='.surface_field.tout'
  integer,parameter 	:: surface_field_output_unit=60
  
  character(len=15)	:: far_field_output_extn='.far_field.fout'
  integer,parameter 	:: far_field_output_unit=61
  
  character(len=30)	:: frequency_output_surface_extn='.frequency_output_surface.fout'
  integer,parameter 	:: frequency_output_surface_unit=62
  
  character(len=9)	:: trcs_output_extn='.rcs.tout'
  integer,parameter 	:: trcs_output_unit=63
  
  character(len=9)	:: rcs_output_extn='.rcs.fout'
  integer,parameter 	:: rcs_output_unit=64
  
  character(len=9)	:: SAR_output_extn='.SAR.fout'
  integer,parameter 	:: SAR_output_unit=65
 
  character(len=10)	:: mode_output_extn='.mode.tout'
  integer,parameter 	:: mode_output_unit=66
  
  character(len=36)	:: frequency_domain_power_surface_extn='.frequency_domain_power_surface.fout'
  integer,parameter 	:: frequency_domain_power_surface_unit=67
  
  character(len=29)	:: frequency_output_volume_extn='.frequency_output_volume.fout'
  integer,parameter 	:: frequency_output_volume_unit=68

  integer,parameter 	:: animation_output_unit=70
  
  character(len=7)	:: record_user_inputs_extn='_in.txt'
  integer,parameter 	:: record_user_inputs_unit=80
  
  integer,parameter 	:: local_file_unit=90
  
  integer,parameter 	:: progress_file_unit=91
  character(len=8)	:: progress_filename='progress'
  
  integer,parameter 	:: scratch_file_unit=92
  
END MODULE File_information
!
! NAME
!     MODULE output_formats
!
! DESCRIPTION
!     declaration of output formats for writing data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/08/2012 CJS
!
!
MODULE output_formats

IMPLICIT NONE

character(len=16),parameter	:: time_domain_output_format='(E16.7,I5,E16.7)'

character(len=20),parameter	:: cable_current_output_format='(E16.7,I5,E16.7,2I5)'

character(len=17),parameter	:: frequency_domain_output_format='(E16.7,I5,5E16.7)'

END MODULE output_formats
!
! NAME
!     MODULE file_header
!
! DESCRIPTION
!     output file header format data
!
!     first lines of data files which identify the format of the following data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 17/08/2012 CJS
!
!
MODULE file_header

IMPLICIT NONE

character(len=18),parameter	:: time_domain_output_file_header='# Time Domain Data'

character(len=23),parameter	:: frequency_domain_output_file_header='# Frequency Domain Data'

END MODULE file_header
