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
! SUBROUTINE calculate_error_Sfilter
! SUBROUTINE calculate_error_Sfilter_PR
! SUBROUTINE calculate_error_Sfilter_PZ
!
! NAME
!     calculate_error_Sfilter
!
! DESCRIPTION
!     calculate the mean square error between filter response (rational function form) and the input data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/1/2013 CJS
!
!
SUBROUTINE calculate_error_Sfilter

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
  
! START

!  CALL write_line('CALLED: calculate_error_Sfilter',0,ff_output_to_screen)

  Mean_square_error=0d0
  
  do freq_loop=1,n_values
  
    do function_loop=1,n_functions

! filter response    
      fit_value=evaluate_Sfilter_frequency_response(filter_S(function_loop),frequency(freq_loop))
      
! add the conductivity response if required
      if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
  
        fit_value=fit_value+filter_sigma(1)/s(freq_loop)
    
      end if   
    
      Mean_square_error=Mean_square_error+abs( fit_value-value(function_loop,freq_loop) )**2
    
    end do ! next function
         
  end do ! next freq_loop
    
  Mean_square_error=Mean_square_error/n_values

!  CALL write_line('FINISHED: calculate_error_Sfilter',0,ff_output_to_screen)

  RETURN
  
END SUBROUTINE calculate_error_Sfilter
!
! NAME
!     calculate_error_Sfilter_PR
!
! DESCRIPTION
!     calculate the mean square error between filter response (pole residue form) and the input data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/1/2013 CJS
!
!
SUBROUTINE calculate_error_Sfilter_PR

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
  
! START

!  CALL write_line('CALLED: calculate_error_Sfilter_PR',0,ff_output_to_screen)

  Mean_square_error=0d0
  
  do freq_loop=1,n_values
  
    do function_loop=1,n_functions

! filter response    
      fit_value=evaluate_Sfilter_PR_frequency_response(filter_S_PR(function_loop),frequency(freq_loop))
      
! add the conductivity response if required
      if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
  
        fit_value=fit_value+filter_sigma(1)/s(freq_loop)
    
      end if   
    
      Mean_square_error=Mean_square_error+abs( fit_value-value(function_loop,freq_loop) )**2
    
    end do ! next function
         
  end do ! next freq_loop
    
  Mean_square_error=Mean_square_error/n_values

!  CALL write_line('FINISHED: calculate_error_Sfilter_PR',0,ff_output_to_screen)

  RETURN
  
END SUBROUTINE calculate_error_Sfilter_PR
!
! NAME
!     calculate_error_Sfilter_PZ
!
! DESCRIPTION
!     calculate the mean square error between filter response (pole zero form) and the input data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/1/2013 CJS
!
!
SUBROUTINE calculate_error_Sfilter_PZ

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
  
! START

!  CALL write_line('CALLED: calculate_error_Sfilter_PZ',0,ff_output_to_screen)

  Mean_square_error=0d0
  
  do freq_loop=1,n_values
  
    do function_loop=1,n_functions

! filter response    
      fit_value=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(function_loop),frequency(freq_loop))
      
! add the conductivity response if required
      if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
  
        fit_value=fit_value+filter_sigma(1)/s(freq_loop)
    
      end if   
    
      Mean_square_error=Mean_square_error+abs( fit_value-value(function_loop,freq_loop) )**2
    
    end do ! next function
         
  end do ! next freq_loop
    
  Mean_square_error=Mean_square_error/n_values

!  CALL write_line('FINISHED: calculate_error_Sfilter_PZ',0,ff_output_to_screen)

  RETURN
  
END SUBROUTINE calculate_error_Sfilter_PZ
