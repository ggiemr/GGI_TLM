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
! SUBROUTINE fourier_transform
!
! NAME
!    fourier_transform
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE fourier_transform

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer		:: function_number

  real*8	:: fmin,fmax,fstep
  real*8	:: log_fmin,log_fmax,log_fstep,log_f
  integer	:: n_frequencies
  integer	:: frequency_loop
  integer	:: timestep
  
  real*8	:: dt
  real*8	:: w
  complex*16	:: integral
  
  character*3	:: freq_range_type
  
! START

  n_functions_of_time=1
  n_functions_of_frequency=1
  
  CALL Allocate_post_data()
  
  function_number=1
  
  CALL read_Time_Domain_Data(function_number) ! read first and only function of time

! setup the frequency list

100 CONTINUE
  write(*,*)'Enter the frequency range type (log or lin)'
  read(*,'(A3)')freq_range_type
  CALL convert_to_lower_case(freq_range_type,3)
  
  if ( (freq_range_type.NE.'log').AND.(freq_range_type.NE.'lin') ) then
    write(*,*)"Frequency range type should be 'log' or 'lin'"
    GOTO 100
  end if
  
  write(record_user_inputs_unit,'(A3)')freq_range_type
  
  write(*,*)'Enter minimum frequency, fmin'
  read(*,*)fmin
  write(*,*)'Enter maximum frequency, fmax'
  read(*,*)fmax
  write(*,*)'Enter the number of frequencies'
  read(*,*)n_frequencies
  
  write(record_user_inputs_unit,'(E16.7,A)')fmin,' fmin'
  write(record_user_inputs_unit,'(E16.7,A)')fmax,' fmax'
  write(record_user_inputs_unit,'(I16,A)')n_frequencies,' n_frequencies'
  
  if (freq_range_type.EQ.'log') then
  
    log_fmin=log10(fmin)
    log_fmax=log10(fmax)
    log_fstep=(log_fmax-log_fmin)/dble(n_frequencies-1)
  
  else if (freq_range_type.EQ.'lin') then
  
    fstep=(fmax-fmin)/dble(n_frequencies-1)
  
  end if

! Generate the frequency list
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
  
    if (freq_range_type.EQ.'log') then
    
      log_f=log_fmin+(frequency_loop-1)*log_fstep

      function_of_frequency(function_number)%frequency(frequency_loop)=10D0**log_f
  
    else if (freq_range_type.EQ.'lin') then
    
      function_of_frequency(function_number)%frequency(frequency_loop)=fmin+(frequency_loop-1)*fstep
  
    end if
   
  end do

! Allocate Fourier transform data  
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  dt=function_of_time(function_number)%time(2)-function_of_time(function_number)%time(1)
  write(*,*)'Timestep=',dt

! calculate the Fourier integral at each of the frequencies specified
  do frequency_loop=1,n_frequencies
    
    w=2d0*pi*function_of_frequency(function_number)%frequency(frequency_loop)
    
    integral=(0d0,0d0)

    do timestep=1,function_of_time(function_number)%n_timesteps
    
      integral=integral+function_of_time(function_number)%value(timestep)	&
                       *exp(-j*w*function_of_time(function_number)%time(timestep))*dt
		       
    end do ! next timestep

    function_of_frequency(function_number)%value(frequency_loop)=integral
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
  
  
END SUBROUTINE fourier_transform
