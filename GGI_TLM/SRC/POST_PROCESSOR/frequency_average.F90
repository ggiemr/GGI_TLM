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
! SUBROUTINE frequency_average
!
! NAME
!    frequency_average
!
! DESCRIPTION
!     
!     
! COMMENTS
!     This should be OK for linear frequency data.
!     We could weight the average so that it would be better for no-linear frequency data
!
! HISTORY
!
!     started 18/01/2013 CJS
!
!
SUBROUTINE frequency_average


USE post_process
USE file_information

IMPLICIT NONE

! local variables

  real*8	:: average
  integer	:: n_average
  
  integer	:: function_number
  integer	:: n_frequencies,n_frequencies_out
  integer	:: frequency_loop
  integer	:: average_loop
  integer	:: average_loop_start
  integer	:: loop
       
  character 	:: bandwidth_type
  character 	:: input_type
  character 	:: average_type
  
  real*8	:: frequency
  real*8	:: percentage_bandwidth
  real*8	:: bandwidth
  real*8	:: fmin_average,fmax_average
  real*8	:: percentage_uniform
  real*8	:: uniform_bandwidth
  real*8	:: fmin_uniform,fmax_uniform
  real*8	:: scaling_factor
  real*8	:: average_factor
  real*8	:: fmin,fmax
  
  real*8 	:: value,magnitude,power_dB
    
! START

  write(*,*)'Average Frequency domain data'
  
  n_functions_of_time=0
  n_functions_of_frequency=2
  
  CALL Allocate_post_data()
  
  write(*,*)'File for frequency averaging'
  function_number=1
  CALL read_frequency_domain_data(function_number)

! read the averaging specification
10 CONTINUE
  write(*,*)"Do you want to specify fixed or percentage bandwidth? ('f' or 'p')"
  read(*,'(A)')bandwidth_type
  CALL convert_to_lower_case(bandwidth_type,1)
  
  if ( (bandwidth_type.NE.'f').AND.(bandwidth_type.NE.'p') ) GOTO 10
  write(record_user_inputs_unit,'(A)')bandwidth_type
  
  if (bandwidth_type.EQ.'f') then
  
    write(*,*)'Enter the averaging bandwidth'
    read(*,*)bandwidth
    write(record_user_inputs_unit,*)bandwidth,' averaging bandwidth'
  
  else ! bandwidth_type.EQ.'p'
  
    write(*,*)'Enter the averaging percentage bandwidth'
    read(*,*)percentage_bandwidth
    write(record_user_inputs_unit,*)percentage_bandwidth,' averaging percentage bandwidth'
  
  end if
  
  write(*,*)'A trapezoidal window is applied in the frequecy average'
  write(*,*)'Enter the percentage of the averaging bandwidth which is uniform'
  read(*,*)percentage_uniform
  write(record_user_inputs_unit,*)percentage_uniform,' uniform window percentage bandwidth'
  
20 CONTINUE
  write(*,*)"Is the input field or power? ('f' or 'p' )"
  read(*,'(A)')input_type
  CALL convert_to_lower_case(input_type,1)
  
  if ( (input_type.NE.'f').AND.(input_type.NE.'p') ) GOTO 20
  write(record_user_inputs_unit,'(A)')input_type

30 CONTINUE
  write(*,*)"Do you want to average field, power or dB? ('f', 'p' or 'd')"
  read(*,'(A)')average_type
  CALL convert_to_lower_case(average_type,1)
  
  if ( (average_type.NE.'f').AND.(average_type.NE.'p').AND.(average_type.NE.'d') ) GOTO 30
  write(record_user_inputs_unit,'(A)')average_type
  
! Do the process twice, the first time to count the numnber of frequency outputs
! and the second to set the frequency output data

  do loop=1,2
  
    n_frequencies=function_of_frequency(1)%n_frequencies
    n_frequencies_out=0
    average_loop_start=1
    
    fmin=function_of_frequency(1)%frequency(1)
    fmax=function_of_frequency(1)%frequency(n_frequencies)
    
    do frequency_loop=1,n_frequencies
  
      frequency=function_of_frequency(1)%frequency(frequency_loop)
      
! work out the frequency range to average over
      if (bandwidth_type.EQ.'p') then  
        bandwidth=frequency*percentage_bandwidth/100d0
      end if
      
      uniform_bandwidth=bandwidth*percentage_uniform/100d0
     
      fmin_average=frequency-bandwidth/2d0
      fmax_average=frequency+bandwidth/2d0

      fmin_uniform=frequency-uniform_bandwidth/2d0
      fmax_uniform=frequency+uniform_bandwidth/2d0
  
      average_factor=0d0
      average=0d0
      n_average=0
      
! check that the whole of the averaging bandwidth sits within the frequency range of the data    
      
      if ( (fmin_average.GE.fmin).AND.(fmax_average.LE.fmax) ) then
      
        do average_loop=average_loop_start,n_frequencies
	
	  if ( (function_of_frequency(1)%frequency(average_loop).GE.fmin_average).AND.	&
	       (function_of_frequency(1)%frequency(average_loop).LE.fmax_average) ) then
! this frequency is in the averaging range

! work out the scaling factor for this sample

            CALL calc_window_function(frequency,fmin_average,fmin_uniform,fmax_uniform,fmax_average,scaling_factor)
            average_factor=average_factor+scaling_factor
	    
            if (n_average.EQ.0) average_loop_start=average_loop ! help to improve efficiency here...
	    
            n_average=n_average+1
	    
	    if (input_type.eq.'f') then
	      value=function_of_frequency(1)%magnitude(average_loop)	    
	    else if (input_type.eq.'p') then
	      value=sqrt( abs(dble(function_of_frequency(1)%value(average_loop))) )
	    end if
	    
	    if (average_type.eq.'f') then ! average type field
	    
	      average=average+value*scaling_factor
	 
	    else if (average_type.eq.'p') then ! average type power
	    
	      average=average+value*value*scaling_factor
	 
	    else if (average_type.eq.'d') then ! average type dB
	    
	      average=average+20d0*log10(value)*scaling_factor
	 
            end if
	    
          end if ! this frequency is in the averaging range
      
          if (function_of_frequency(1)%frequency(average_loop).GT.fmax_average) GOTO 1000 ! finish loop
      
        end do ! next frequency
	
1000    CONTINUE ! jump here if the frequency goes above fmax_average
      
      end if ! whole of the averaging bandwidth sits within the frequency range of the data 
      
      if (n_average.ne.0) then
      
! we can calculate an average
!!!        average=average/(n_average*average_factor)
        average=average/average_factor
        n_frequencies_out=n_frequencies_out+1
	
        if (loop.eq.2) then
	
	  if (average_type.eq.'f') then 

	    magnitude=average
	    power_dB=20d0*log10(magnitude)
	 
	  else if (average_type.eq.'p') then 
	 
	    magnitude=sqrt(average)
	    power_dB=20d0*log10(magnitude)
	 
	  else if (average_type.eq.'d') then 
	  
	    magnitude=10d0**(average/20d0)
	    power_dB=average
	 
	  end if
	  	
	  if (input_type.eq.'f') then
	  
   	    function_of_frequency(2)%frequency(n_frequencies_out)	=frequency
   	    function_of_frequency(2)%value(n_frequencies_out)	=(0d0,0d0) ! only average mag and dB
   	    function_of_frequency(2)%magnitude(n_frequencies_out) =magnitude
   	    function_of_frequency(2)%phase(n_frequencies_out)	=0d0       ! only average mag and dB
   	    function_of_frequency(2)%dB(n_frequencies_out)	=power_dB
	    	  
	  else if (input_type.eq.'p') then
	  
   	  function_of_frequency(2)%frequency(n_frequencies_out)	=frequency
   	  function_of_frequency(2)%value(n_frequencies_out)	=cmplx(magnitude**2)
   	  function_of_frequency(2)%magnitude(n_frequencies_out) =magnitude**2
   	  function_of_frequency(2)%phase(n_frequencies_out)	=0d0       ! only average mag and dB
   	  function_of_frequency(2)%dB(n_frequencies_out)	=power_dB
	    
	  end if
	 
        end if ! loop=2 so set data
	
      end if ! we have an average at this frequency
        
    end do ! next frequency
  
    if (loop.eq.1) then
    
      write(*,*)'Number of average frequencies=',n_frequencies_out
    
      function_of_frequency(2)%n_frequencies=n_frequencies_out
      ALLOCATE ( function_of_frequency(2)%frequency(1:n_frequencies_out) )
      ALLOCATE ( function_of_frequency(2)%value(1:n_frequencies_out) )
      ALLOCATE ( function_of_frequency(2)%magnitude(1:n_frequencies_out) )
      ALLOCATE ( function_of_frequency(2)%phase(1:n_frequencies_out) )
      ALLOCATE ( function_of_frequency(2)%dB(1:n_frequencies_out) )
      
    end if
  
  end do ! next loop

  function_number=2
  CALL write_Frequency_Domain_Data(function_number)

  CALL Deallocate_post_data()

  RETURN
  

  
END SUBROUTINE frequency_average
!
! __________________________________________________
!
!
SUBROUTINE calc_window_function(f,f1,f2,f3,f4,value)

real*8 f,f1,f2,f3,f4,value

! START

  if ( (f.GT.f1).AND.(f.LE.f2) ) then
     
    value=(f-f1)/(f2-f1)
    
  else if ( (f.GT.f2).AND.(f.LE.f3) ) then
     
    value=1d0  
    
  else if ( (f.GT.f3).AND.(f.LE.f4) ) then
        
    value=(f-f4)/(f3-f4)
  
  else
  
    value=0d0  
    
  end if
  
END SUBROUTINE calc_window_function

