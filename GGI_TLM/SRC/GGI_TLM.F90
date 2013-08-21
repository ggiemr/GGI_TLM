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
! PROGRAM GGI_TLM
!
! NAME
!     GGI_TLM
!
! DESCRIPTION
!     TLM solver
!     The solver runs on the already created meshed model from MODEL_BUILDER and MODEL_CHECKER
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!     Parallel implementation started 21/11/2012 CJS
!
!
PROGRAM GGI_TLM

USE TLM_general
USE geometry

IMPLICIT NONE

! local variables

  logical :: read_data_for_computation_only

! START

#if defined(SEQ)

  np=1
  rank=0
  
#elif defined(MPI)

! initialise MPI
  CALL MPI_INIT(ierr)
  
! get the rank of this processor
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)  
  
! get the total number of processors (np)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
  
#endif
  
  write(*,*)'Rank',rank,' number of processes= ',np

  CALL write_progress('STARTED: GGI_TLM')
  
  CALL write_line('GGI_TLM',0,output_to_screen_flag)
    
  CALL write_license()

  CALL reset_packet_data()
  
  CALL read_problem_name()
  
  CALL open_general_files()
  
  CALL reset_TLM_solver_parameters()
  
  CALL read_mesh()
  
!  CALL trim_mesh() ! redundant with new mesh generation...
  
  read_data_for_computation_only=.TRUE.
  CALL read_cable_model(read_data_for_computation_only)
  
  CALL read_input_file_solver()
  
  CALL check_solver_input_data()
  
  CALL set_TLM_solver_parameters()

  CALL allocate_temporary_mesh_arrays()
  
  CALL set_volume_material_mesh()
  
  CALL set_surface_material_mesh()
 
  CALL set_outputs_in_mesh()
  
  CALL set_cable_bundles_in_mesh()
  
  if (np.gt.1) CALL partition_cable_mesh()
  
  CALL initialise_mode_stir_surfaces()
  
  CALL set_huygens_surface_data()
  
  CALL set_excitations_in_mesh() ! flags the excitations in local_cell_excitation or local_surface_excitation
  
  CALL set_cell_update_codes()
  
  CALL set_face_update_codes()
  
  CALL initialise_outputs()
  
  CALL initialise_cable_bundles()
  
  CALL initialise_excitations()

  CALL deallocate_temporary_mesh_arrays()
  
  CALL initialise_excitation_functions()

  CALL allocate_mesh()
  
  CALL calculate_volume_material_filter_coefficients
  
  CALL calculate_surface_material_filter_coefficients
  
  CALL allocate_volume_material_filter_data()
  
  CALL allocate_surface_material_filter_data()
  
  CALL initialise_wrap_outer_boundary()

  CALL run_TLM()
  
  CALL write_frequency_domain_outputs()

  CALL finish_outputs()

  CALL deallocate_mesh()

  CALL deallocate_geometry()
  
  CALL deallocate_materials()
  
  CALL deallocate_outputs()
  
  CALL deallocate_excitations()
  
  CALL deallocate_mode_stir()
  
  CALL deallocate_cables()
  
  CALL write_line('FINISHED: GGI_TLM',0,output_to_screen_flag)
  
  CALL write_line_integer('Number of warnings=',number_of_warnings,0,output_to_screen_flag)
    
#if defined(MPI)

  CALL MPI_FINALIZE(ierr)
  
#endif

  CALL write_progress('FINISHED: GGI_TLM')
  
END PROGRAM GGI_TLM
