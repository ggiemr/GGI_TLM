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
! SUBROUTINE RMS_frequency_domain_data
!
! NAME
!    RMS_frequency_domain_data
!
! DESCRIPTION
!     Calculate the Square root of the Sum of squares of complex data or alternatively the RMS value
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/11/2014 CJS based on RMS_frequency_domain_d
!
!
SUBROUTINE RMS_frequency_domain_data


USE post_process
USE file_information

IMPLICIT NONE

! local variables
  
  complex*16	:: f1,f2,result
  
  integer	:: function_number
  integer	:: n_functions
  integer	:: n_frequencies
  integer	:: frequency_loop
  
  character	:: ch
  logical	:: RMS_flag
  
! START
 
  write(*,*)'Do you want to calculate the square root of the sum of squares (s) or RMS value (r) ?'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  if (ch.eq.'s') then
    RMS_flag=.FALSE.
  else if (ch.eq.'r') then
    RMS_flag=.TRUE.
  else
    write(*,*)"ERROR: response should be 's' or 'r'"
    STOP   
  end if
 
  write(*,*)' '
  write(*,*)'Result=sqrt( sum 1..N {fn*conjg(fn)) )'
  write(*,*)' '
  write(*,*)'where fn is a set of N frequency domain quantities'
  write(*,*)' '
  
  write(*,*)'Enter the number of frequency domain quantities to sum:'
  read(*,*)n_functions
  write(record_user_inputs_unit,*)n_functions,' number of functions to sum'
  
  write(post_process_info_unit,*)'	Number of functions to sum=',n_functions

  n_functions_of_time=0
  n_functions_of_frequency=n_functions+1
  
  CALL Allocate_post_data()
  
  write(post_process_info_unit,*)'	Frequency domain functions:'
  do function_number=1,n_functions
  
    write(*,*)'File for frequency domain data:'
    CALL read_frequency_domain_data(function_number)
     
    if (function_number.gt.1) then
! check that the frequencies match...
      if (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(function_number)%n_frequencies) then
        write(*,*)'Frequency mismatch in functions f 1 and f',function_number
        write(*,*)'n_frequencies, f1=',function_of_frequency(1)%n_frequencies
        write(*,*)'n_frequencies, fn=',function_of_frequency(function_number)%n_frequencies
        STOP
      end if
  
      do frequency_loop=1,function_of_frequency(1)%n_frequencies
  
        if ( function_of_frequency(1)%frequency(frequency_loop).NE.	&
             function_of_frequency(function_number)%frequency(frequency_loop) ) then
	 
          write(*,*)'Frequency mismatch in functions f 1 and f',function_number
          write(*,*)'frequency number',frequency_loop
          write(*,*)'frequency, f1=',function_of_frequency(1)%frequency(frequency_loop)
          write(*,*)'frequency, fn=',function_of_frequency(function_number)%frequency(frequency_loop)
          STOP
      
        end if
  
      end do
      
    end if !function_number.gt.1
    
  end do ! next function number

! Allocate memory for result

  n_frequencies=function_of_frequency(1)%n_frequencies
  function_number=n_functions_of_frequency
  
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                function_of_frequency(1)%frequency(1:n_frequencies)
     
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  
  do frequency_loop=1,n_frequencies
  
    result=(0d0,0d0)
    
    do function_number=1,n_functions

      result=result+( function_of_frequency(function_number)%value(frequency_loop)*	&
              conjg(function_of_frequency(function_number)%value(frequency_loop)) )
      
    end do
    
    if (RMS_flag) then
      result=sqrt(result/cmplx(n_functions))
    else
      result=sqrt(result)
    end if
    
    function_number=n_functions_of_frequency
    
    function_of_frequency(function_number)%value(frequency_loop)=result
    
    function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(function_number)%value(frequency_loop))
    function_of_frequency(function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(function_number)%value(frequency_loop))   )
    function_of_frequency(function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))

  end do ! next frequency value
  
  function_number=n_functions_of_frequency
  CALL write_Frequency_Domain_Data(function_number)

  CALL Deallocate_post_data()

  RETURN
  

  
END SUBROUTINE RMS_frequency_domain_data
