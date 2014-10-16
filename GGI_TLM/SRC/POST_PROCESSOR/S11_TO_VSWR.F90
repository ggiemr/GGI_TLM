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
! SUBROUTINE S11_TO_VSWR
!
! NAME
!    S11_TO_VSWR
!
! DESCRIPTION
!     calculate VSWR from complex S11 data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/10/2014 CJS  for horn antenna test case comparison with
!                            manufacturer's data
!
SUBROUTINE S11_TO_VSWR


USE post_process
USE file_information

IMPLICIT NONE

! local variables
  
  real*8	:: f
  complex*16	:: S11
  real*8	:: rho,VSWR
 
  integer	:: n_frequencies
  integer	:: frequency_loop
  integer	:: function_number
    
! START

!  write(*,*)'S11_TO_VSWR transformation'

  n_functions_of_time=0
  n_functions_of_frequency=2
  
  CALL Allocate_post_data()
  
  write(*,*)'Read S11 data'
  CALL read_frequency_domain_data(1)
    
  n_frequencies=function_of_frequency(1)%n_frequencies

  function_number=2
    
  function_of_frequency(function_number)%n_frequencies=n_frequencies  
    
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		    &
  		function_of_frequency(1)%frequency(1:n_frequencies)	
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
    
    S11=function_of_frequency(1)%value(frequency_loop)
    rho=abs(S11)
    VSWR= (1d0+rho)/(1d0-rho)
        
    function_of_frequency(2)%value(frequency_loop)=VSWR
    function_of_frequency(2)%magnitude(frequency_loop)= 	  &
        	  abs(function_of_frequency(2)%value(frequency_loop))
    function_of_frequency(2)%phase(frequency_loop)=   &
        	    atan2( imag(function_of_frequency(2)%value(frequency_loop)), &
        		   dble(function_of_frequency(2)%value(frequency_loop))       )
    function_of_frequency(2)%dB(frequency_loop)=      &
        	    20d0*log10(function_of_frequency(2)%magnitude(frequency_loop))

  end do ! next frequency value
  
  write(*,*)'Write VSWR data to file'
  CALL write_Frequency_Domain_Data(2)

  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE S11_TO_VSWR

