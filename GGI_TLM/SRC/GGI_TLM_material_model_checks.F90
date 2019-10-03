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
! PROGRAM GGI_TLM_material_model_checks
!
! NAME
!     GGI_TLM_material_model_checks
!
! DESCRIPTION
!     TLM material model checks:
!     Allows the visualisation of meshes on a material type basis
!     Calculation of material frequency response
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!     17/10/2014 CJS - add write_all_material_info_to_file flag and options 5 and 6:
!                      eliminates the need for user input if information on all materials is required
!                      This option is used in the automatic documentation
!
!
PROGRAM GGI_TLM_material_model_checks

USE TLM_general
USE geometry
USE TLM_surface_materials
USE TLM_volume_materials

IMPLICIT NONE

! local variables

  integer :: number_of_options,option
  
  logical :: write_all_material_info_to_file

! START

  CALL write_progress('STARTED: GGI_TLM_material_model_checks')
  
  CALL write_line('GGI_TLM_material_model_checks',0,output_to_screen_flag)
  
  CALL write_license()
  
  CALL reset_packet_data()
  
  CALL read_problem_name()

  CALL open_general_files()
  
  CALL reset_TLM_solver_parameters()
  
  CALL read_mesh()
  
!  CALL trim_mesh() ! redundant with new mesh generation...
  
  CALL read_input_file_solver()
  
  CALL check_solver_input_data()

  CALL allocate_temporary_mesh_arrays()
  
  CALL set_volume_material_mesh()
  
  CALL set_pml_volume_material_mesh()
  
  CALL set_surface_material_mesh()

  
  write(*,*)'Number of volume materials=',n_volume_materials
  write(*,*)'Number of volumes=',n_volumes
  write(*,*)'Number of surface materials=',n_surface_materials
  write(*,*)'Number of surfaces=',n_surfaces
  
  number_of_options=6

10 CONTINUE  ! start of post_processing action

  write(*,*)
  write(*,*)'Material model checking options are:'
  write(*,*)
  write(*,*)'1. View volume material frequency response '
  write(*,*)'2. View volume material cells '
  write(*,*)'3. View surface material frequency response '
  write(*,*)'4. View surface material faces '
  write(*,*)'5. Write all volume material information to file '
  write(*,*)'6. Write all surface material information to file '
  write(*,*)
  
  write(*,'(A,I2,A)')'Please enter the required cable model option 1 :',number_of_options,' or 0 to quit'
  read(*,*)option
  
  if (option.EQ.0) then  ! close files, deallocate memory and stop

    CALL deallocate_temporary_mesh_arrays()
  
    CALL deallocate_mesh()

    CALL deallocate_geometry()
  
    CALL deallocate_materials()

    CALL write_progress('FINISHED: GGI_TLM_material_model_checks')
  
    CALL write_line('FINISHED: GGI_TLM_material_model_checks',0,output_to_screen_flag)

    STOP
     
  else if (option.EQ.1) then
    
    write(*,*)'View volume material frequency response'
    
    write_all_material_info_to_file=.FALSE.
    CALL Volume_material_frequency_response(write_all_material_info_to_file)
    
  else if (option.EQ.2) then
    
    write(*,*)'View volume material cells'
    
    write_all_material_info_to_file=.FALSE.
    CALL plot_volume_material_cells(write_all_material_info_to_file)      
    CALL plot_mesh_PML_volumes()
    CALL plot_mesh_boundary()

  else if (option.EQ.3) then
  
    write(*,*)'View surface material frequency response'
       
    write_all_material_info_to_file=.FALSE.
    CALL Surface_material_frequency_response(write_all_material_info_to_file)

  else if (option.EQ.4) then
  
    write(*,*)'View surface material faces'
    
    write_all_material_info_to_file=.FALSE.
    CALL plot_surface_material_faces(write_all_material_info_to_file)
    CALL plot_mesh_boundary()

  else if (option.EQ.5) then
  
    write(*,*)'Write all volume material information to file'
    
    write_all_material_info_to_file=.TRUE.
    CALL Volume_material_frequency_response(write_all_material_info_to_file)
    CALL plot_volume_material_cells(write_all_material_info_to_file)  
    CALL plot_mesh_PML_volumes()
    CALL plot_mesh_boundary()

  else if (option.EQ.6) then
  
    write(*,*)'Write all surface material information to file'
    
    write_all_material_info_to_file=.TRUE.
    CALL Surface_material_frequency_response(write_all_material_info_to_file)
    CALL plot_surface_material_faces(write_all_material_info_to_file)
    CALL plot_mesh_boundary()

  end if
  
  GOTO 10
    
  
END PROGRAM GGI_TLM_material_model_checks
