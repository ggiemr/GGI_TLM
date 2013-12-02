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
! SUBROUTINE stabilise_filter
! SUBROUTINE stabilise_dielectric_material
! SUBROUTINE stabilise_magnetic_material
! SUBROUTINE stabilise_thin_layer
! SUBROUTINE stabilise_impedance
!
! NAME
!     stabilise_filter
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
SUBROUTINE stabilise_filter

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE filter_functions

IMPLICIT NONE
 
! local variables

  integer	:: function_loop
  integer	:: i
  
! START

! 1. Calculate the mean square error of the initial filter function  

  CALL calculate_error_Sfilter()  
  CALL write_line_real('Initial Sfilter.             Mean square error=',Mean_square_error,0,ff_output_to_screen)

! 2. Transform the filter into Pole-Residue format  

  do function_loop=1,n_functions
    filter_S_PR(function_loop)=Convert_filter_S_to_S_PR( filter_S(function_loop) )
  end do ! next function
  
! 3. Calculate the mean square error of the Pole-Residue format filter as a check on the transformation
 
  CALL calculate_error_Sfilter_PR()
  CALL write_line_real('Initial Sfilter_PR.          Mean square error=',Mean_square_error,0,ff_output_to_screen)

! 4. Stabilise the poles of the filter functions i.e. Re(pole)<0
                  
  do function_loop=1,n_functions
    do i=1,filter_S_PR(function_loop)%order
      if (dble(filter_S_PR(function_loop)%poles(i)).gt.0d0) then 
!  	write(*,*)'Stabilise pole',i,function_loop
  	filter_S_PR(function_loop)%poles(i)=filter_S_PR(function_loop)%poles(i)-2d0*dble(filter_S_PR(function_loop)%poles(i))
  	filter_S_PR(function_loop)%residues(i)=-filter_S_PR(function_loop)%residues(i)  ! ***** also change sign of residue... *****
      end if
    end do
  end do 

! 5. Calculate the mean square error of the stabilised Pole-Residue format filter 
 
  CALL calculate_error_Sfilter_PR()
  CALL write_line_real('Pole Stabilised Sfilter_PR,  Mean square error=',Mean_square_error,0,ff_output_to_screen)

!6. Call the appropriate stabilisation routine for the fit-type

  if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
    
    CALL stabilise_dielectric_or_magnetic_material()
       
  else if (fit_type.eq.thin_layer) then  
    
    CALL stabilise_thin_layer()
    
  else if (fit_type.eq.impedance) then  
    
    CALL stabilise_impedance()
    
  else if (fit_type.eq.general) then  
    
! No action
    
  end if

! 7. Calculate the mean square error of the fully stabilised Pole-Residue format filter 
 
  CALL calculate_error_Sfilter_PR()
  CALL write_line_real('Final Stabilised Sfilter_PR, Mean square error=',Mean_square_error,0,ff_output_to_screen)

! 8. Transform the filter back into Rational function format  

  do function_loop=1,n_functions
    filter_S(function_loop)=Convert_filter_S_PR_to_S( filter_S_PR(function_loop) )
  end do ! next function
  
! 9. Calculate the mean square error of the Pole-Residue format filter as a check on the transformation
 
  CALL calculate_error_Sfilter()
  CALL write_line_real('Final Stabilised Sfilter.    Mean square error=',Mean_square_error,0,ff_output_to_screen)


END SUBROUTINE stabilise_filter
