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
! SUBROUTINE read_input_file_solver
!
! NAME
!     read_input_file_solver
!
! DESCRIPTION
!     read input file packets related to the TLM solver
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
SUBROUTINE read_input_file_solver

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

character*256	:: input_line

! START
    
  rewind(input_file_unit)
  
10  CONTINUE
    
! read line from input file
    read(input_file_unit,'(A)',end=100)input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,256)

    if (input_line.EQ.'simulation_time') then
    
      CALL read_simulation_time()
     
    else if (input_line.EQ.'volume_material_list') then
    
      CALL read_volume_material_list()
     
    else if (input_line.EQ.'surface_material_list') then
    
      CALL read_surface_material_list()
      
    else if (input_line.EQ.'outer_boundary_reflection_coefficient') then
    
      CALL read_outer_boundary_reflection_coefficient()    
     
    else if (input_line.EQ.'excitation_function_list') then
    
      CALL read_excitation_function_list()
     
    else if (input_line.EQ.'excitation_point_list') then
    
      CALL read_excitation_point_list()
     
    else if (input_line.EQ.'excitation_surface_list') then
    
      CALL read_excitation_surface_list()
     
    else if (input_line.EQ.'huygens_surface') then
    
      CALL read_huygens_surface()
     
    else if (input_line.EQ.'output_point_list') then
    
      CALL read_output_point_list()
     
    else if (input_line.EQ.'output_surface_list') then
    
      CALL read_output_surface_list()
     
    else if (input_line.EQ.'output_volume_list') then
    
      CALL read_output_volume_list()
     
    else if (input_line.EQ.'output_volume_average_list') then
    
      CALL read_output_volume_average_list()
     
    else if (input_line.EQ.'sar_volume_list') then
    
      CALL read_sar_volume_list()
     
    else if (input_line.EQ.'far_field_surface') then
    
      CALL read_far_field_surface()
     
    else if (input_line.EQ.'frequency_output_surface_list') then
    
      CALL read_frequency_output_surface()
     
    else if (input_line.EQ.'frequency_output_volume_list') then
    
      CALL read_frequency_output_volume()
     
    else if (input_line.EQ.'frequency_domain_power_surface_list') then
                            
      CALL read_frequency_domain_power_surface()
     
    else if (input_line.EQ.'rcs_surface') then
    
      CALL read_rcs_surface()
     
    else if (input_line.EQ.'excitation_mode_list') then
    
      CALL read_excitation_mode_list()
     
    else if (input_line.EQ.'output_mode_list') then
    
      CALL read_output_mode_list()
            
    else if (input_line.EQ.'mode_stir_surface_list') then
    
      CALL read_mode_stir_surface_list()
            
    else if (input_line.EQ.'bicubic_warp_flag') then
    
      CALL read_bicubic_warp_flag()
            
    else if (input_line.EQ.'frequency_scale_flag') then
    
      CALL read_frequency_scale_flag()
            
    else if (input_line.EQ.'wrapping_boundary_conditions') then
    
      CALL read_wrapping_boundary_conditions()
      
    end if
    
    GOTO 10
    
100 CONTINUE
  
  CALL write_line('FINISHED: read_input_file_solver',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE read_input_file_solver
