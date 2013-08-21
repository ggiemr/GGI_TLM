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
! SUBROUTINE get_testing_frequencies
!
! NAME
!     get_testing_frequencies
!
! DESCRIPTION
!     
! Get frequencies for stability testing:
! Use the input frequencies then add some more frequencies at higher and lower frequencies 
! in order to make the stability testing more robust
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
SUBROUTINE get_testing_frequencies

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_testing_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE
 
! local variables

  real*8	:: fmax,fmin
  
  integer	:: freq_scale_loop
  integer	:: freq_count
  integer	:: freq_loop

! START

  CALL write_line('CALLED: get_testing_frequencies',0,ff_output_to_screen)

  n_testing_frequencies=16+n_values+16
  
  ALLOCATE( testing_frequency(1:n_testing_frequencies) )

  fmin=frequency(1)
  fmax=frequency(n_values)
  
  freq_count=0
  
! add some low frequency terms on a logarithmic scale

  do freq_scale_loop=16,1,-1
      
    freq_count=freq_count+1
    testing_frequency(freq_count)=fmin/(2d0**freq_scale_loop)
	
  end do

! include all the input frequencies

  do freq_loop=1,n_values

    freq_count=freq_count+1
    testing_frequency(freq_count)=frequency(freq_loop)
   
  end do ! next freq_loop

  
! add some high frequency terms on a logarithmic scale

  do freq_scale_loop=1,16
      
    freq_count=freq_count+1
    testing_frequency(freq_count)=fmax*(2d0**freq_scale_loop)
	
  end do
  
! construct the normalised complex testing frequency array  
  ALLOCATE (  testing_s(1:n_testing_frequencies)  )
  testing_s(1:n_testing_frequencies)=2d0*pi*j*testing_frequency(1:n_testing_frequencies)/wnorm
  
! write the testing frequencies to file
 
  OPEN(unit=testing_frequency_file_unit,file=testing_frequency_filename)
  
  do freq_count=1,n_testing_frequencies
    write(testing_frequency_file_unit,8000)freq_count,testing_frequency(freq_count)
8000 format(I10,E16.6)
  end do
  
  CLOSE(unit=testing_frequency_file_unit)


  CALL write_line('FINISHED: get_testing_frequencies',0,ff_output_to_screen)

END SUBROUTINE get_testing_frequencies
