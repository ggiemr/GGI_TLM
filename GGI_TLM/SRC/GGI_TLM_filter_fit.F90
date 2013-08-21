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
!
! NAME
!     PROGRAM GGI_TLM_filter_fit
!
! DESCRIPTION
!     filter function fitting process based on the weiner hopf method
!     plus stabilisation then optimisation for developing 
!     frequency dependent material/ thin layer/ impedance models for GGI_TLM
!     and other codes which represent frequency dependent parameters as
!     rational functions in complex frequency. 
!
!     Process:
!	1. Read initial data
!	2. Stabilise initial data to give the input data to the model fitting process (and write to file)
!       3. Create initial filter function(s) using the Weiner Hopf method with no stability constraints
!       4. Stabilise the filter function(s)
!       5. Optimise the filter function(s) whilst maintaining stability
!       6. Write filter function(s) to file together with the filter frequency response(s)
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 11/12/12 CJS
!     based on process developed in the Flaviir project and HIRF-SE at the 
!     University of Nottingham
!
!
  PROGRAM GGI_TLM_filter_fit
!
       
USE TLM_general
USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE filter_types
USE filter_functions

IMPLICIT NONE
       
! local variables

  integer	:: function_loop

! START

  CALL write_progress('STARTED: GGI_TLM_filter_fit')

  ff_output_to_screen=.TRUE.
  
  CALL write_line('GGI_TLM_filter_fit',0,output_to_screen_flag)
  
  CALL write_license()
       
! read problem name, data and options associated with this filter fit      
  CALL read_filter_fit_information()

! sort out any stability issues with the input data      
  CALL stabilise_FF_input_data()
       
! Write the stabilised input data to file. This is the data to which the model will be fitted
  CALL write_FF_input_data()

! get a set of frequencies for testing the stability of the filter model.  
  CALL get_testing_frequencies()
  
! get an initial set of filter coefficients by the Weiner Hopf method.
  CALL get_initial_filter_coefficients()
  
! test the stability of this filter
  CALL test_stability()
  
  if (stabilise_filter_flag) then
! Ensure that this filter is stable
    CALL stabilise_filter()
  end if
  
! test the stability of the stabilised filter
  CALL test_stability()
  
  if (optimise_filter_flag) then
! optimise the filter coefficients
    CALL optimise_filter()
  end if
  
! test the stability of the optimised filter
  CALL test_stability()
  
! write the filter coefficients to file
  CALL write_filter()
  
! write the filter frequency response to file
  CALL write_filter_frequency_response()
  
  CALL deallocate_filter_fit_memory()

  CALL write_progress('FINISHED: GGI_TLM_filter_fit')

  END
       
