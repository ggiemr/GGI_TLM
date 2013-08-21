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
! SUBROUTINE combine_frequency_domain_data
! SUBROUTINE combine_frequency_domain_magnitude_data
!
! NAME
!    combine_frequency_domain_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 18/01/2013 CJS
!
!
SUBROUTINE combine_frequency_domain_data


USE post_process
USE file_information

IMPLICIT NONE

! local variables

  real*8	:: A,B,C
  
  complex*16	:: f1,f2,result
  
  integer	:: function_number
  integer	:: n_frequencies
  integer	:: frequency_loop
  
! START

!  write(*,*)'Combine Frequency domain data'
  
  write(*,*)' '
  write(*,*)'Result=((A*f1)/(B*f2)+C'
  write(*,*)' '
  write(*,*)'where f1 and f2 are frequency domain quantities and A, B and C are constants'
  write(*,*)' '

  n_functions_of_time=0
  n_functions_of_frequency=3
  
  CALL Allocate_post_data()
  
  write(*,*)'File for f1 data:'
  function_number=1
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'File for f2 data:'
  function_number=2
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Enter constant A'
  read(*,*)A
  write(*,*)'Enter constant B'
  read(*,*)B
  write(*,*)'Enter constant C'
  read(*,*)C
  
  write(record_user_inputs_unit,*)A,' constant, A'
  write(record_user_inputs_unit,*)B,' constant, B'
  write(record_user_inputs_unit,*)C,' constant, C'

! check that the frequencies match...

  if (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(2)%n_frequencies) then
    write(*,*)'Frequency mismatch in functions f1 and f2'
    write(*,*)'n_frequencies, f1=',function_of_frequency(1)%n_frequencies
    write(*,*)'n_frequencies, f2=',function_of_frequency(2)%n_frequencies
    STOP
  end if
  
  do frequency_loop=1,function_of_frequency(1)%n_frequencies
  
    if ( function_of_frequency(1)%frequency(frequency_loop).NE.	&
         function_of_frequency(2)%frequency(frequency_loop) ) then
	 
      write(*,*)'Frequency mismatch in functions f1 and f2'
      write(*,*)'frequency number',frequency_loop
      write(*,*)'frequency, f1=',function_of_frequency(1)%frequency(frequency_loop)
      write(*,*)'frequency, f2=',function_of_frequency(2)%frequency(frequency_loop)
      STOP
      
    end if
  
  end do

! Allocate memory for result

  n_frequencies=function_of_frequency(1)%n_frequencies
  function_number=3
  
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                function_of_frequency(1)%frequency(1:n_frequencies)
     
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
  
    f1=function_of_frequency(1)%value(frequency_loop)
    f2=function_of_frequency(2)%value(frequency_loop)
    
    result=(A*f1)/(B*f2)+C
    
    function_of_frequency(function_number)%value(frequency_loop)=result
    
    function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(function_number)%value(frequency_loop))
    function_of_frequency(function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(function_number)%value(frequency_loop))   )
    function_of_frequency(function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))

  end do ! next frequency value
  
  CALL write_Frequency_Domain_Data(function_number)

  CALL Deallocate_post_data()

  RETURN
  

  
END SUBROUTINE combine_frequency_domain_data
!
! NAME
!    combine_frequency_domain_magnitude_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 21/01/2013 CJS
!
!
SUBROUTINE combine_frequency_domain_magnitude_data


USE post_process
USE file_information

IMPLICIT NONE

! local variables

  real*8	:: A,B,C
  
  real*8	:: f1,f2,result
  
  integer	:: function_number
  integer	:: n_frequencies
  integer	:: frequency_loop
  
! START

!  write(*,*)'Combine Frequency domain magnitude data'
  
  write(*,*)' '
  write(*,*)'Result=((A*|f1|)/(B*|f2|)+C'
  write(*,*)' '
  write(*,*)'where f1 and f2 are frequency domain quantities and A, B and C are constants'
  write(*,*)' '

  n_functions_of_time=0
  n_functions_of_frequency=3
  
  CALL Allocate_post_data()
  
  write(*,*)'File for f1 data:'
  function_number=1
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'File for f2 data:'
  function_number=2
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Enter constant A'
  read(*,*)A
  write(*,*)'Enter constant B'
  read(*,*)B
  write(*,*)'Enter constant C'
  read(*,*)C
  
  write(record_user_inputs_unit,*)A,' constant, A'
  write(record_user_inputs_unit,*)B,' constant, B'
  write(record_user_inputs_unit,*)C,' constant, C'

! check that the frequencies match...

  if (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(2)%n_frequencies) then
    write(*,*)'Frequency mismatch in functions f1 and f2'
    write(*,*)'n_frequencies, f1=',function_of_frequency(1)%n_frequencies
    write(*,*)'n_frequencies, f2=',function_of_frequency(2)%n_frequencies
    STOP
  end if
  
  do frequency_loop=1,function_of_frequency(1)%n_frequencies
  
    if ( function_of_frequency(1)%frequency(frequency_loop).NE.	&
         function_of_frequency(2)%frequency(frequency_loop) ) then
	 
      write(*,*)'Frequency mismatch in functions f1 and f2'
      write(*,*)'frequency number',frequency_loop
      write(*,*)'frequency, f1=',function_of_frequency(1)%frequency(frequency_loop)
      write(*,*)'frequency, f2=',function_of_frequency(2)%frequency(frequency_loop)
      STOP
      
    end if
  
  end do

! Allocate memory for result

  n_frequencies=function_of_frequency(1)%n_frequencies
  function_number=3
  
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                function_of_frequency(1)%frequency(1:n_frequencies)
     
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
  
    f1=function_of_frequency(1)%magnitude(frequency_loop)
    f2=function_of_frequency(2)%magnitude(frequency_loop)
    
    result=(A*f1)/(B*f2)+C
    
    function_of_frequency(function_number)%value(frequency_loop)=dcmplx(result,0d0)
    
    function_of_frequency(function_number)%magnitude(frequency_loop)=result
    function_of_frequency(function_number)%phase(frequency_loop)=0d0
    function_of_frequency(function_number)%dB(frequency_loop)=20d0*log10(result)

  end do ! next frequency value
  
  CALL write_Frequency_Domain_Data(function_number)

  CALL Deallocate_post_data()

  RETURN
  

  
END SUBROUTINE combine_frequency_domain_magnitude_data
