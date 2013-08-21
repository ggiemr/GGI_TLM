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
! SUBROUTINE read_input_file_cables
!
! NAME
!     read input file cables
!
! DESCRIPTION
!     read cable information from the input file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
SUBROUTINE read_input_file_cables

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

character*256	:: input_line

! START
  
  CALL write_line('CALLED: read_input_file_cables',0,output_to_screen_flag)
    
  rewind(input_file_unit)
  
10  CONTINUE
    
! read line from input file
    read(input_file_unit,'(A)',end=100)input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,256)

    if (input_line.EQ.'cable_geometry_list') then
    
      CALL read_cable_geometry_list()
          
    else if (input_line.EQ.'cable_list') then
    
      CALL read_cable_list()
          
    else if (input_line.EQ.'cable_junction_list') then
    
      CALL read_cable_junction_list()
          
    else if (input_line.EQ.'cable_output_list') then
    
      CALL read_cable_output_list()
    
    else if (input_line.EQ.'surface_material_list') then
    
      CALL read_surface_material_list()  ! this is required so as to indicate whether a surface 
                                         ! has material properties set or not.

! The following packets redefine solution parameters and are dealt with here      

    else if (input_line.EQ.'number_of_Fourier_terms_in_pul_lc_calc') then
      
      read(input_file_unit,*)number_of_Fourier_terms_in_PUL_LC_calc

    else if (input_line.EQ.'tlm_cell_equivalent_radius_factor') then
      
      read(input_file_unit,*)TLM_cell_equivalent_radius_factor

    else if (input_line.EQ.'capacitance_equivalent_radius_factor') then
      
      read(input_file_unit,*)Capacitance_equivalent_radius_factor

    else if (input_line.EQ.'inductance_equivalent_radius_factor') then
      
      read(input_file_unit,*)Inductance_equivalent_radius_factor

    else if (input_line.EQ.'max_cable_bundle_diameter_factor') then
      
      read(input_file_unit,*)Max_cable_bundle_diameter_factor

    else if (input_line.EQ.'lc_correction_type_geometry_scale') then
      
      Cable_LC_Correction_type=LC_correction_type_geometry_scale

    else if (input_line.EQ.'lc_correction_type_subtract_cell_inductance') then
      
      Cable_LC_Correction_type=LC_correction_type_subtract_cell_inductance
     
    end if
    
    GOTO 10
    
100 CONTINUE
  
  CALL write_line('FINISHED: read_input_file_cables',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE read_input_file_cables
