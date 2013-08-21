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
! SUBROUTINE write_filter_frequency_response
!
! NAME
!     write_filter_frequency_response
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/12/2012 CJS
!
!
SUBROUTINE write_filter_frequency_response

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
  integer	:: freq_loop
  
  complex*16    :: fit_value
  
  character(len=256) :: temp_filename
  character(len=256) :: filename(1:4)

! START

  CALL write_line('CALLED: write_filter_frequency_response',0,ff_output_to_screen)
  
  CALL add_integer_to_filename(FF_name,order,temp_filename)

  if (fit_type.eq.dielectric_material) then 
  
    filename(1)=trim(temp_filename)//dielectric_trial_extension
 
  else if (fit_type.eq.magnetic_material) then 
  
    filename(1)=trim(temp_filename)//magnetic_trial_extension
 
  else if (fit_type.eq.thin_layer) then 
  
    filename(1)=trim(temp_filename)//thin_layer_z11_trial_extension
    filename(2)=trim(temp_filename)//thin_layer_z12_trial_extension
    filename(3)=trim(temp_filename)//thin_layer_z21_trial_extension
    filename(4)=trim(temp_filename)//thin_layer_z22_trial_extension
 
  else if (fit_type.eq.impedance) then 
   
    filename(1)=trim(temp_filename)//impedance_trial_extension

  end if
  
  do function_loop=1,n_functions
  
    OPEN(unit=filter_fit_data_file_unit,file=filename(function_loop))
  
    do freq_loop=1,n_values

! filter response    
      fit_value=evaluate_Sfilter_frequency_response(filter_S(function_loop),frequency(freq_loop))
      
! add the conductivity response if required
      if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
  
        fit_value=fit_value+filter_sigma(1)/s(freq_loop)
    
      end if   
  
      write(filter_fit_data_file_unit,8000)frequency(freq_loop),real(fit_value),imag(fit_value)
8000  format(3E16.6)
     
    end do ! next freq_loop
    
  end do ! next function

  CALL write_line('FINISHED: write_filter_frequency_response',0,ff_output_to_screen)

  RETURN
  
END SUBROUTINE write_filter_frequency_response
