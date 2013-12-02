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
!
! SUBROUTINE FFT_MAIN
!
! NAME
!    FFT
!
! DESCRIPTION
!     controlling subroutine for the FFT
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/11/2013 CJS
!
!
SUBROUTINE FFT_MAIN

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer		:: function_number

  real*8	:: fmin,fmax,fstep
  integer	:: n_frequencies
  integer	:: frequency_loop
  integer	:: timestep,n_timesteps
  
  integer 	:: i
  
  real*8	:: dt
  
  complex*16,allocatable	:: x(:)
  
  character 	:: ch  
  logical	:: dt_scale_flag
  
! START

  n_functions_of_time=1
  n_functions_of_frequency=1
  
  CALL Allocate_post_data()
  
  function_number=1
  
  CALL read_Time_Domain_Data(function_number) ! read first and only function of time

! setup the frequency limits

  n_timesteps=function_of_time(function_number)%n_timesteps
  dt=function_of_time(function_number)%time(2)-function_of_time(function_number)%time(1)
  
5 CONTINUE       
  write(*,*)'Do you want to scale the FFT by the timestep? (y or n)'
  read(*,'(A)')ch
  if ( (ch.eq.'y').OR.(ch.eq.'Y') ) then
    dt_scale_flag=.TRUE.
  else if ( (ch.eq.'n').OR.(ch.eq.'N') ) then
    dt_scale_flag=.FALSE.
  else
    write(*,*)"Response should be 'y' or 'n'"
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')ch

! set the number of frequencies to be equal to the next power of 2 above the number of timesteps  

  i=1
10  CONTINUE
    i=i*2
    if (i.ge.n_timesteps) then
      n_frequencies=i
    else
      GOTO 10
    end if

  fmin=0d0
  fstep=1d0/(n_frequencies*dt)
  fmax=fstep*(n_frequencies-1)

! Generate the frequency list
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
      
    function_of_frequency(function_number)%frequency(frequency_loop)=fmin+(frequency_loop-1)*fstep
  
  end do

! Allocate Fourier transform data  
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  write(*,*)'Timestep=',dt
  
! copy the time domain data into the temporary complex array, x

  ALLOCATE( x(1:n_frequencies) )
  
  x(1:n_timesteps)=function_of_time(function_number)%value(1:n_timesteps)
  
  if (n_timesteps.LT.n_frequencies) then
    x(n_timesteps+1:n_frequencies)=(0d0,0d0)
  end if
  
  CALL FFT(x,n_frequencies)
  
! Put the FFT data into the appropriate format
  do frequency_loop=1,n_frequencies
  
    if (dt_scale_flag) then
      function_of_frequency(function_number)%value(frequency_loop)=dt*x(frequency_loop)
    else
      function_of_frequency(function_number)%value(frequency_loop)=x(frequency_loop)
    end if
    
    function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(function_number)%value(frequency_loop))
    function_of_frequency(function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(function_number)%value(frequency_loop))   )
    function_of_frequency(function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))

  end do ! next frequency value
  
  CALL write_Frequency_Domain_Data(function_number)

  DEALLOCATE( x )
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE FFT_MAIN
!
! SUBROUTINE FFT
!
! NAME
!    FFT
!
! DESCRIPTION
!     Fast Fourier Transform routine - called recursively
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
RECURSIVE SUBROUTINE FFT(x,N)

USE constants

IMPLICIT NONE

  integer	:: n
  complex*16	:: x(n)

! local variables

  integer			:: n2
  complex*16,allocatable	:: xe(:)
  complex*16,allocatable	:: xo(:)
  complex*16	:: const
  
  integer	:: i

! START

! No action required for n=1  
  if(n .LE. 1) RETURN
  
  n2=n/2
 
  ALLOCATE(xe(1:n2))
  ALLOCATE(xo(1:n2))
 
! fill odd and even data
 
  do i=1,n2
    xo(i)=x(2*i-1)
    xe(i)=x(2*i)
  end do

! FFT odd and even sequences
  CALL FFT(xo,n2)
  CALL FFT(xe,n2)
 
! combine odd and even FFTs

  do i=1,n2
 
    const=exp(-2d0*pi*j*(i-1)/n)

    x(i)   = xo(i)+xe(i)*const
    x(i+N2)= xo(i)-xe(i)*const
    
  end do
 
  DEALLOCATE(xo)
  DEALLOCATE(xe)
  
  RETURN
 
END SUBROUTINE FFT
