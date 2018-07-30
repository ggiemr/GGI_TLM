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
! SUBROUTINE Apply_filter_in_frequency_domain
! SUBROUTINE Calculate_filter_frequency_response
!
! NAME
!    Apply_filter_in_frequency_domain
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/1/2013 CJS
!
!
SUBROUTINE Apply_filter_in_frequency_domain

USE post_process
USE filter_types		        
USE filter_functions		        
USE file_information

IMPLICIT NONE

! local variables

  integer			:: function_number

! Filter stuff
  character*256			:: filter_file_name
  type(Sfilter) 		:: Sfilter1
  type(Zfilter) 		:: Zfilter1

  real*8			:: frequency
  complex*16			:: f1,f2,result
  
  integer	:: n_frequencies
  integer	:: frequency_loop
  character 	:: ch
  logical	:: Zfilter_flag
  real*8	:: dt

! START

! Allocate and read the time domain input data
  n_functions_of_time=0
  n_functions_of_frequency=2
  
  CALL Allocate_post_data()
  
  function_number=1
      
  CALL Read_Frequency_Domain_Data(function_number) ! read function of frequency 
  
! Allocate and read the filter function to apply

  write(*,*)'Enter the filename for the filter function'
  read(*,'(A)')filter_file_name

  open(UNIT=local_file_unit, 					    &
       FILE=trim(filter_file_name),	    &
       STATUS='old',							    &
       ERR=9000)

  call read_Sfilter(Sfilter1,local_file_unit) ! filter function

  close(UNIT=local_file_unit)
  
  write(record_user_inputs_unit,'(A)')trim(filter_file_name)  
  write(post_process_info_unit,*)'	Filter filename:',trim(filter_file_name)

! we offer a choice of applying the filter frequency response as 
! specified or that of the digital filter obtained using the bilinear transformation  
  write(*,*)"Enter 'S' to apply the S domain filter frequency response as specified "
  write(*,*)"or 'Z' to apply the digital filter response obtained using the "
  write(*,*)"bilinear transformation for a given timestep"
 
  read(*,'(A)')ch
  if ( (ch.eq.'s').OR.(ch.eq.'S') ) then
    Zfilter_flag=.FALSE.
    write(post_process_info_unit,*)'	Apply s-plane filter'
  else if ( (ch.eq.'z').OR.(ch.eq.'Z') ) then
    Zfilter_flag=.TRUE.
    write(post_process_info_unit,*)'	Apply z-plane filter'
  else
    write(*,*)"Response should be 'u' or 'n'"
    STOP
  end if
  write(record_user_inputs_unit,'(A)')ch
  
  if (Zfilter_flag) then
    write(*,*)'Enter the timestep'
    read(*,*)dt
    write(record_user_inputs_unit,*)dt,' timestep'
    write(post_process_info_unit,*)'	Timestep=',dt,' seconds'
    Zfilter1=s_to_z(Sfilter1,dt) 
  end if
  
! Apply the S or Z domain filter to the frequency domain data

  n_frequencies=function_of_frequency(1)%n_frequencies
  
  function_number=2
! Allocate memory for result
  
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                function_of_frequency(1)%frequency(1:n_frequencies)
     
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  function_of_frequency(function_number)%n_frequencies=n_frequencies
  
  do frequency_loop=1,n_frequencies
  
    frequency=function_of_frequency(1)%frequency(frequency_loop)    
    function_of_frequency(function_number)%frequency(frequency_loop)=frequency
    
    f1=function_of_frequency(1)%value(frequency_loop)   
     
    if (Zfilter_flag) then
      f2=evaluate_Zfilter_frequency_response(Zfilter1,frequency)    
    else    
      f2=evaluate_Sfilter_frequency_response(Sfilter1,frequency)    
    end if
    
    result=f1*f2
        
    function_of_frequency(function_number)%value(frequency_loop)=result
    
    function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(function_number)%value(frequency_loop))
    function_of_frequency(function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(function_number)%value(frequency_loop))   )
    function_of_frequency(function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))

  end do ! next frequency value 
  
! Write the filtered data set to file
  
  CALL write_frequency_domain_data(function_number)

! deallocate memory  
  CALL Deallocate_post_data()
  CALL deallocate_Sfilter( Sfilter1 )
  
  RETURN
  
9000 CALL write_line('Error reading filter function',0,.TRUE.)
     CALL write_line('Problem opening filter file:'//trim(filter_file_name),0,.TRUE.)
     STOP
  
END SUBROUTINE Apply_filter_in_frequency_domain
!
! NAME
!    Calculate_filter_frequency_response
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/1/2013 CJS
!
!
SUBROUTINE Calculate_filter_frequency_response

USE post_process
USE filter_types		        
USE filter_functions		        
USE file_information

IMPLICIT NONE

! local variables

  integer			:: function_number

! Filter stuff
  character*256			:: filter_file_name
  type(Sfilter) 		:: Sfilter1
  type(Zfilter) 		:: Zfilter1

  real*8			:: frequency
  complex*16			:: result
  
  integer	:: n_frequencies
  integer	:: frequency_loop
    
  character*3	:: freq_range_type
  real*8	:: fmin,fmax,fstep
  real*8	:: log_fmin,log_fmax,log_fstep,log_f
  
  character 	:: ch
  logical	:: Zfilter_flag
  real*8	:: dt
  
  logical,parameter :: warp_flag=.FALSE.
  real*8,parameter  :: warp_scale=1d0

! START

! Allocate and read the time domain input data
  n_functions_of_time=0
  n_functions_of_frequency=1
  
  CALL Allocate_post_data()
  
! Allocate and read the filter function to apply

  write(*,*)'The filter specification should be in the form of the following example:'
  write(*,*)'1.0E+08  # wnormalisation constant' 
  write(*,*)'1        # a order' 
  write(*,*)'4.0 0.75 # a coefficients' 
  write(*,*)'1        # b order' 
  write(*,*)'1.0 0.5  # b coefficients' 
  write(*,*)'' 

  write(*,*)'Enter the name of file containing the filter specification'
  read(*,'(A)')filter_file_name
  write(record_user_inputs_unit,'(A)')trim(filter_file_name)

  open(UNIT=local_file_unit, 					    &
       FILE=trim(filter_file_name),	    &
       STATUS='old',							    &
       ERR=9000)

  call read_Sfilter(Sfilter1,local_file_unit) ! filter function

  close(UNIT=local_file_unit)

! we offer a choice to calculate the filter frequency response as 
! specified or that of the digital filter obtained using the bilinear transformation  
  write(*,*)"Enter 'S' to apply the S domain filter frequency response as specified "
  write(*,*)"or 'Z' to apply the digital filter response obtained using the "
  write(*,*)"bilinear transformation for a given timestep"
 
  read(*,'(A)')ch
  if ( (ch.eq.'s').OR.(ch.eq.'S') ) then
    Zfilter_flag=.FALSE.
  else if ( (ch.eq.'z').OR.(ch.eq.'Z') ) then
    Zfilter_flag=.TRUE.
  else
    write(*,*)"Response should be 's' or 'z'"
    STOP
  end if
  write(record_user_inputs_unit,'(A)')ch
  
  if (Zfilter_flag) then
    write(*,*)'Enter the timestep'
    read(*,*)dt
    write(record_user_inputs_unit,*)dt,' timestep'
    Zfilter1=s_to_z_warp(Sfilter1,dt,warp_flag,warp_scale) 
  end if
  
! Apply the S or Z domain filter to the frequency domain data

100 CONTINUE

  function_number=1

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
  
! Evaluate the S domain filter response in the frequency domain
  
  function_number=1
  
  do frequency_loop=1,n_frequencies
  
    frequency=function_of_frequency(1)%frequency(frequency_loop)    
    function_of_frequency(function_number)%frequency(frequency_loop)=frequency
      
    if (Zfilter_flag) then
      result=evaluate_Zfilter_frequency_response(Zfilter1,frequency)    
    else    
      result=evaluate_Sfilter_frequency_response(Sfilter1,frequency)    
    end if
        
    function_of_frequency(function_number)%value(frequency_loop)=result
    
    function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(function_number)%value(frequency_loop))
    function_of_frequency(function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(function_number)%value(frequency_loop))   )
    function_of_frequency(function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))

  end do ! next frequency value 
  
! Write the output data set to file
  
  CALL write_frequency_domain_data(function_number)

! deallocate memory  
  CALL Deallocate_post_data()
  CALL deallocate_Sfilter( Sfilter1 )
  
  RETURN
  
9000 CALL write_line('Error reading filter function',0,.TRUE.)
     CALL write_line('Problem opening filter file:'//trim(filter_file_name),0,.TRUE.)
     STOP
  
END SUBROUTINE Calculate_filter_frequency_response
!
! NAME
!    Calculate_filter_time_response
!
! DESCRIPTION
!    Take a filter function in the S plane, apply the bilinear transformation (with timestep or TLM cell size specified)
!    The calculate the impulse response of the resulting Z plane filter
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 23/9/2014 CJS
!
!
SUBROUTINE Calculate_filter_time_response

USE post_process
USE filter_types		        
USE filter_functions		        
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer			:: function_number

! Filter stuff
  character*256			:: filter_file_name
  type(Sfilter) 		:: Sfilter1
  type(Zfilter) 		:: Zfilter1
  type(Zfilter_response) 	:: Zfilter_data1

  integer			:: timestep,n_timesteps
  real*8			:: f_input

  real*8	:: dl,dt
  
  character	:: ch
  
  logical,parameter :: warp_flag=.FALSE.
  real*8,parameter  :: warp_scale=1d0

! START

! Allocate and read the time domain input data
  n_functions_of_time=1
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
! Allocate and read the filter function to apply

  write(*,*)'The filter specification should be in the form of the following example:'
  write(*,*)'1.0E+08  # wnormalisation constant' 
  write(*,*)'1        # a order' 
  write(*,*)'4.0 0.75 # a coefficients' 
  write(*,*)'1        # b order' 
  write(*,*)'1.0 0.5  # b coefficients' 
  write(*,*)'' 

  write(*,*)'Enter the name of file containing the filter specification'
  read(*,'(A)')filter_file_name
  write(record_user_inputs_unit,'(A)')trim(filter_file_name)

  open(UNIT=local_file_unit, 					    &
       FILE=trim(filter_file_name),	    &
       STATUS='old',							    &
       ERR=9000)

  call read_Sfilter(Sfilter1,local_file_unit) ! filter function

  close(UNIT=local_file_unit)

! we offer a choice to calculate the filter frequency response as 
! specified or that of the digital filter obtained using the bilinear transformation  
  write(*,*)"Enter 't' to set the timestep or "
  write(*,*)"or 'l' to specify a TLM cell size and calculate the timestep appropriate for this"
 
  read(*,'(A)')ch
  
  write(record_user_inputs_unit,'(A)')ch
  
  if ( (ch.eq.'t').OR.(ch.eq.'T') ) then
  
    write(*,*)'Enter the timestep'
    read(*,*)dt
    write(record_user_inputs_unit,*)dt,' timestep' 
    
  else if ( (ch.eq.'l').OR.(ch.eq.'L') ) then
  
    write(*,*)'Enter the TLM cell size'
    read(*,*)dl
    write(record_user_inputs_unit,*)dl,' TLM cell size' 
    
    dt=dl/(2d0*c0)

  else
    write(*,*)"Response should be 't' or 'l'"
    STOP
  end if

! Calculate the Z domain filter
  Zfilter1=s_to_z_warp(Sfilter1,dt,warp_flag,warp_scale) 

  Zfilter_data1=allocate_Zfilter_response(Zfilter1%a%order,Zfilter1%b%order)

! allocate memory for the time domain response function

  write(*,*)'Enter the number of timesteps of the impulse response to calculate'
  read(*,*)n_timesteps
  write(record_user_inputs_unit,*)n_timesteps,' n_timesteps' 
  
  function_of_time(1)%n_timesteps=n_timesteps
  ALLOCATE ( function_of_time(1)%time(1:n_timesteps) )
  ALLOCATE ( function_of_time(1)%value(1:n_timesteps) )
  
! Apply the digital filter to the time domain data
  
  Zfilter_data1%w(:)=0d0
  
  do timestep=1,n_timesteps
  
    CALL timeshift_Zfilter(Zfilter_data1)
    
! create impulse driving function
    f_input=0d0
    if (timestep.eq.1) then
      f_input=1d0
    end if
    
    CALL evaluate_Zfilter(Zfilter1,Zfilter_data1,f_input)
    
    function_of_time(1)%time(timestep)=dble(timestep-1)*dt
    function_of_time(1)%value(timestep)=Zfilter_data1%f  
    
  end do ! next sample
  
! Write the filtered data set to file
  function_number=1
  
  CALL write_time_domain_data(function_number)

! deallocate memory  
  CALL Deallocate_post_data()
  CALL deallocate_Sfilter( Sfilter1 )
  CALL deallocate_Zfilter( Zfilter1 )
  CALL deallocate_Zfilter_data( Zfilter_data1 )
  
  RETURN
  
9000 CALL write_line('Error reading filter function',0,.TRUE.)
     CALL write_line('Problem opening filter file:'//trim(filter_file_name),0,.TRUE.)
     STOP
  
END SUBROUTINE Calculate_filter_time_response
