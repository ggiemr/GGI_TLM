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
!     stabilise_impedance
!
! DESCRIPTION
!     Given the impedance filter function in pole-residue form:
!     1. Ensure that all poles are stable (LHS of s plane)
!     2. Ensure that Re{Z}>=0 for all testing frequencies
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/12/2012 CJS
!
!
SUBROUTINE stabilise_impedance

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_testing_data
USE FF_file_stuff
USE filter_functions

IMPLICIT NONE
 
! local variables
       
  real*8	:: min_real_Z
 
  integer	:: i,term
  real*8	:: Re_fs
  complex*16 	:: fs_fit

! START
      	
! check whether the constant term is positive i.e. the limit as f-> infinity
  if(filter_S_PR(1)%C.lt.0d0) then   
!    print*,'Stabilising filter_S_PR(1)%C:',filter_S_PR(1)%C
    filter_S_PR(1)%C=-filter_S_PR(1)%C/10d0
  end if
       
  min_real_Z=0d0
		  
! loop over frequencies in testing list	
  do i=1,n_testing_frequencies

! evaluate Z(jw) at this frequency	

    fs_fit=evaluate_Sfilter_PR_frequency_response(filter_S_PR(1),testing_frequency(i))
	  
! test that real part is greater than zero	  
    Re_fs=(dble(fs_fit))

    if (Re_fs.lt.min_real_Z) min_real_Z=Re_fs
	 
  end do ! next testing frequency

  filter_S_PR(1)%C=filter_S_PR(1)%C-min_real_Z*1.1D0  ! add to constant term to ensure that the impedance is positive

  RETURN

END SUBROUTINE stabilise_impedance
