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
! PROGRAM GGI_TLM_cable_model_checks
!
! NAME
!     GGI_TLM_cable_model_checks
!
! DESCRIPTION
!     Check cable model before running the TLM solver:
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/08/2012 CJS
!
!
PROGRAM GGI_TLM_cable_model_checks

USE TLM_general
USE geometry
USE Cables

IMPLICIT NONE

! local variables

  integer :: number_of_options,option

  logical :: read_data_for_computation_only

! START
  
  CALL write_progress('STARTED: GGI_TLM_cable_model_checks')
  
  CALL write_line('GGI_TLM_cable_model_checks',0,output_to_screen_flag)
  
  CALL write_license()
  
  CALL read_problem_name()
  
  CALL read_mesh()
  
!  CALL trim_mesh() ! redundant with new mesh generation...

  read_data_for_computation_only=.FALSE.
  CALL read_cable_model(read_data_for_computation_only)
  
  write(*,*)'Number of cable geometries=',n_cable_geometries
  write(*,*)'Number of cables=',n_cables
  write(*,*)'Number of bundle segments=',n_bundle_segments
  write(*,*)'Number of bundle segment geometries=',n_bundle_segment_geometries
  write(*,*)'Number of cell centre junctions=',n_cell_centre_junctions 
  write(*,*)'Number of face junctions=',n_face_junctions
  write(*,*)'Number of cable outputs=',n_cable_outputs
  
  number_of_options=5

10 CONTINUE  ! start of post_processing action

  write(*,*)
  write(*,*)'Cable model checking options are:'
  write(*,*)
  write(*,*)'1. Output the cable geometry  '
  write(*,*)'2. Output the cable route (segment based) '
  write(*,*)'3. Output the cable bundle segment geometry at any segment '
  write(*,*)'4. Output the cable junction specification at any cell '
  write(*,*)'5. Output the cable face junction specification at any cell face '
  write(*,*)
  
  write(*,'(A,I2,A)')'Please enter the required cable model option 1 :',number_of_options,' or 0 to quit'
  read(*,*)option
  
  if (option.EQ.0) then  ! close files, deallocate memory and stop

    CALL deallocate_mesh()

    CALL deallocate_geometry()
  
    CALL deallocate_materials()
  
    CALL deallocate_outputs()
  
    CALL deallocate_excitations()
  
    CALL deallocate_cables()
  
    CALL write_progress('FINISHED: GGI_TLM_cable_model_checks')
  
    CALL write_line('FINISHED: GGI_TLM_cable_model_checks',0,output_to_screen_flag)
  
    STOP
    
  else if (option.EQ.1) then
    
    write(*,*)'Output the cable geometry '
    
    CALL Output_cable_geometry()
    
  else if (option.EQ.2) then
    
    write(*,*)'Output the cable route (segment based) '
    
    CALL Output_cable_route()
    
  else if (option.EQ.3) then
  
    write(*,*)'Output the cable bundle geometry at any segment '
    
    CALL Output_bundle_geometry()

  else if (option.EQ.4) then
  
    write(*,*)'Output the cable junction specification at any cell '
    
    CALL Output_junction_specification()
    
  else if (option.EQ.5) then
  
    write(*,*)'Output the cable face junction specification at any cell face '
    
    CALL Output_face_junction_specification()
 
  end if
  
  GOTO 10
 
END PROGRAM GGI_TLM_cable_model_checks
