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
! PROGRAM GGI_TLM_cable_model_builder
!
! NAME
!     GGI_TLM_cable_model_builder
!
! DESCRIPTION
!     GGI_TLM_cable_model_builder solver
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
PROGRAM GGI_TLM_cable_model_builder

USE TLM_general
USE geometry
USE File_information

IMPLICIT NONE

! local variables

! START

  CALL write_line('GGI_TLM_cable_model_builder',0,output_to_screen_flag)
  
  CALL write_license()
  
  CALL write_progress('STARTED: GGI_TLM_cable_model_builder')
  
  CALL reset_packet_data()
  
  CALL read_problem_name()

  CALL open_general_files()

  CALL open_file(cable_info_file_unit,cable_info_file_extn)
  
  CALL read_mesh()
   
!  CALL trim_mesh() ! redundant with new mesh generation

  CALL read_input_file_cables()

  CALL check_cable_input_data()
  
  CALL set_TLM_solver_parameters() ! this sets the timestep
				    
  CALL create_cable_LCRG_matrices()  ! calculate internal L,C for shielded cables 
  
  CALL create_cable_meshes()  ! set the cable to mesh segment list: cable_list(cable)%cable_segment_list(segment)
				
  CALL create_bundle_segments() ! bundle_segment_list(bundle_segment_count)%cable_list(1:number_of_cables)

  CALL get_bundle_segment_geometries() ! work out the number of different bundle segment geometries
				    
  CALL create_bundle_LCRG_matrices()   ! calculate the L, C and R matrices for the different bundle segment geometries

  CALL build_cell_centre_junction_list() ! allocate and start to populate the cell_centre_junction_list
                                         ! on a cable by cable basis
  CALL build_face_junction_list() 

  CALL set_bundle_excitations() 

  CALL set_bundle_outputs() 

  CALL write_cable_model()     

  CALL deallocate_geometry()
  
  CALL deallocate_cables()

  CALL close_file(cable_info_file_unit)
  
  CALL write_progress('FINISHED: GGI_TLM_cable_model_builder')
  
  CALL write_line('FINISHED: GGI_TLM_cable_model_builder',0,output_to_screen_flag)
  
END PROGRAM GGI_TLM_cable_model_builder
