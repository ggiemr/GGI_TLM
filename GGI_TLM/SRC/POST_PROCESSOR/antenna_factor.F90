! SUBROUTINE complex_antenna_factor
! SUBROUTINE complex_antenna_factor_2
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
!    complex_antenna_factor
!
! DESCRIPTION
!     Calculate complex antenna factor from S parameter (S21) measurement data for two identical antennas
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/10/2013 CJS
!
!
SUBROUTINE complex_antenna_factor


USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  real*8	:: d
  real*8	:: Zload
  real*8	:: lambda
  
  real*8	:: mag,phase,last_phase
  
  integer	:: n_twopi
  
  complex*16	:: S21,beta,prop,AFsqr,AF
  
  integer	:: function_number
  integer	:: n_frequencies
  integer	:: frequency_loop
  
! START

!  write(*,*)'Complex_antenna_factor'
  
  n_functions_of_time=0
  n_functions_of_frequency=2
  
  CALL Allocate_post_data()
  
  write(*,*)'File for complex S21 data:'
  write(post_process_info_unit,*)'	Complex S21 data:'
  function_number=1
  CALL read_frequency_domain_data(function_number)
    
  write(*,*)'Enter antenna separation (m)'
  read(*,*)d
  
  write(record_user_inputs_unit,*)d,' Antenna separation, d'
  write(post_process_info_unit,*)'	Antenna separation=',d,' m'
    
  write(*,*)'Enter antenna load impedance (ohms)'
  read(*,*)Zload
  
  write(record_user_inputs_unit,*)Zload,' Antenna load impedance, Zload'
  write(post_process_info_unit,*)'	Antenna load impedance=',Zload,' ohms'

! Allocate memory for result

  n_frequencies=function_of_frequency(1)%n_frequencies
  function_number=2
  
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                function_of_frequency(1)%frequency(1:n_frequencies)
     
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
  
    S21=function_of_frequency(1)%value(frequency_loop)
    
    lambda=c0/function_of_frequency(1)%frequency(frequency_loop)
    
    beta=2d0*pi/lambda
    
    prop=(exp(-j*beta*d))/d
    
    AFsqr=(j*z0/(Zload*lambda))*prop/S21
    
! square root AFsqr attempting to get the correct root...

    mag=abs(AFsqr)
    
    phase=atan2( imag(AFsqr),dble(AFsqr) )

! Attempt to make the phase continuous
    if (frequency_loop.gt.1) then
      if ( abs(last_phase-phase).gt.pi) then
        n_twopi=NINT((last_phase-phase)/(2d0*pi))
        phase=phase+2d0*pi*n_twopi
      end if
      
      if ( (last_phase-phase).gt.pi) then
        write(*,*)'Error1 in phase calculation',phase,last_phase
      else if ( (phase-last_phase).gt.pi) then
        write(*,*)'Error2 in phase calculation',phase,last_phase
      end if
      
    end if
   
    last_phase=phase
    
    mag=sqrt(mag)
    phase=phase/2d0
    
    AF=mag*cos(phase)+j*mag*sin(phase)
        
    function_of_frequency(function_number)%value(frequency_loop)=AF
    
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

END SUBROUTINE complex_antenna_factor
!
! NAME
!    complex_antenna_factor_2
!
! DESCRIPTION
!     Calculate complex antenna factor from S parameter (S21) measurement data with one known and one unknown antenna
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/10/2013 CJS
!
!
SUBROUTINE complex_antenna_factor_2


USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  real*8	:: d
  real*8	:: Zload
  real*8	:: lambda
  
  real*8	:: mag,phase,last_phase
  
  real*8	:: df
  
  integer	:: n_twopi
  
  complex*16	:: S21,beta,prop,AF,AF_product,AF_known
  
  integer	:: function_number
  integer	:: n_frequencies
  integer	:: frequency_loop
  
! START

!  write(*,*)'Complex_antenna_factor_2'
  
  n_functions_of_time=0
  n_functions_of_frequency=3
  
  CALL Allocate_post_data()
  
  write(*,*)'File for complex S21 data:'
  write(post_process_info_unit,*)'	Complex S21 data:'
  function_number=1
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'File for complex antenna factor data for known antenna:'
  write(post_process_info_unit,*)'	Complex antenna factor of known antenna:'
  function_number=2
  CALL read_frequency_domain_data(function_number)
    
  write(*,*)'Enter antenna separation (m)'
  read(*,*)d
  
  write(record_user_inputs_unit,*)d,' Antenna separation, d'
  write(post_process_info_unit,*)'	Antenna separation=',d,' m'
    
  write(*,*)'Enter antenna load impedance (ohms)'
  read(*,*)Zload
  
  write(record_user_inputs_unit,*)Zload,' Antenna load impedance, Zload'
  write(post_process_info_unit,*)'	Antenna load impedance=',Zload,' ohms'

! check that the frequencies match...

  if (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(2)%n_frequencies) then
    write(*,*)'Frequency mismatch in functions f1 and f2'
    write(*,*)'n_frequencies, f1=',function_of_frequency(1)%n_frequencies
    write(*,*)'n_frequencies, f2=',function_of_frequency(2)%n_frequencies
    STOP
  end if
  
  do frequency_loop=1,function_of_frequency(1)%n_frequencies
  
    df=function_of_frequency(1)%frequency(frequency_loop)*1D-6
    
    if ( abs(function_of_frequency(1)%frequency(frequency_loop)-	&
             function_of_frequency(2)%frequency(frequency_loop)).GT.df ) then
	 
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
  
    S21=function_of_frequency(1)%value(frequency_loop)
    
    lambda=c0/function_of_frequency(1)%frequency(frequency_loop)
    
    beta=2d0*pi/lambda
    
    prop=(exp(-j*beta*d))/d
    
    AF_product=(j*z0/(Zload*lambda))*prop/S21
    
! The Antenna factor of the unknown antenna is the product of the antenna factors 
! divided by the antenna factor of the known antenna

    AF_known=function_of_frequency(2)%value(frequency_loop)

    AF=AF_product/AF_known
       
    function_of_frequency(function_number)%value(frequency_loop)=AF
    
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

END SUBROUTINE complex_antenna_factor_2
