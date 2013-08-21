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
! SUBROUTINE IELF
!
! NAME
!    IELF
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
SUBROUTINE IELF


USE post_process
USE file_information

IMPLICIT NONE

! local variables

  real*8	:: IELF_value
  
  real*8	:: fmin,fmax
  
  integer	:: function_number
  integer	:: n_samples1,n_samples2
  
  integer	:: frequency_loop
  
  character 	:: average_type
  
  character(len=256)	:: filename
  
  integer	:: sample_min,sample_max,last_sample,sample
  
  real*8	:: freq_1,freq_2
  real*8	:: value_1,value_2,frequency
        
  real*8 	:: power_dB1,power1,field1
  real*8 	:: power_dB2,power2,field2
  
  real*8 	:: error
 
  real*8	:: local_ielf_value
  
  character*20	:: judgement
  
  integer	:: i
  
! START

  write(*,*)'IELF Analysis'
  
  write(*,*)' '

  n_functions_of_time=0
  n_functions_of_frequency=2
  
  CALL Allocate_post_data()
  
  write(*,*)'File 1 with reference frequency data'
  function_number=1
  CALL read_frequency_domain_data(function_number)  
  n_samples1=function_of_frequency(1)%n_frequencies
  
  write(*,*)'File 2 with frequency data to be interpolated'
  function_number=2
  CALL read_frequency_domain_data(function_number)
  n_samples2=function_of_frequency(2)%n_frequencies
  
20 CONTINUE
  write(*,*)"Do you want to process field, power or power(db)? ('f', 'p' or 'd')"
  read(*,'(A)')average_type
  CALL convert_to_lower_case(average_type,1)
  
  if ( (average_type.NE.'f').AND.(average_type.NE.'p').AND.(average_type.NE.'d') ) GOTO 20
  write(record_user_inputs_unit,'(A)')average_type
     
  write(*,*)'Enter the minimum frequency for IELF calculation'
  read(*,*)fmin
  write(record_user_inputs_unit,*)fmin,' minimum frequency for IELF calculation'
  
  write(*,*)'Enter the maximum frequency for IELF calculation'
  read(*,*)fmax
  write(record_user_inputs_unit,*)fmax,' maximum frequency for IELF calculation'

! loop over samples of file 1 to work out sample_min and sample_max for local_ielf_value

  sample_min=0
  sample_max=0
  
  do sample=2,n_samples1-1
  
    if ( ( function_of_frequency(1)%frequency(sample-1).ge.fmin ).AND.( sample_min.eq.0 ) ) sample_min=sample
    if ( ( function_of_frequency(1)%frequency(sample+1).ge.fmax ).AND.( sample_max.eq.0 ) ) sample_max=sample
    
  end do
  
  if (sample_min.eq.0) sample_min=2
  if (sample_max.eq.0) sample_max=n_samples1-1
  
  write(*,*)'n_samples1 =',n_samples1
  write(*,*)'Sample_min=',sample_min,' fmin=',function_of_frequency(1)%frequency(sample_min)
  write(*,*)'Sample_max=',sample_max,' fmax=',function_of_frequency(1)%frequency(sample_max)
       
! loop over samples in frequency range and do local_ielf_value calculation

  local_ielf_value=0d0

  last_sample=1

  do sample=sample_min,sample_max
  
    value_1=function_of_frequency(1)%magnitude(sample)
    frequency=function_of_frequency(1)%frequency(sample)
    
! find the frequencies which lie either side of f1 and interpolate to give the value from file 2
    do i=last_sample,n_samples2-1
    
	 freq_1=function_of_frequency(2)%frequency(i)
	 freq_2=function_of_frequency(2)%frequency(i+1)
	 
      if ( (freq_1.le.frequency).AND.	&
           (freq_2.gt.frequency)  ) then
	   
    	value_2=function_of_frequency(2)%magnitude(i)+						&
	   ( (frequency-freq_1)/(freq_2-freq_1) )*					&
	     (function_of_frequency(2)%magnitude(i+1)-function_of_frequency(2)%magnitude(i))
	    
    	last_sample=i
    	
    	GOTO 1000
	
      end if
      
    end do  ! next sample of file 2

    write(*,*)'Sample not found in file 2, frequency=',frequency
    write(*,*)'First frequency=',function_of_frequency(2)%frequency(last_sample)
    write(*,*)'Last frequency=',function_of_frequency(2)%frequency(n_samples2-1)

    STOP

1000  CONTINUE

! we now have value_1 and value_2 at the current frequency 

! scale values appropriately here...
! input type field
    
    field1=value_1
    power1=field1*field1
    power_dB1=10d0*log10(power1)
    field2=value_2
    power2=field2*field2
    power_dB2=10d0*log10(power2)
   
    if (average_type.eq.'f') then ! average type field

      value_1=field1
      value_2=field2

    else if (average_type.eq.'p') then ! average type power

      value_1=power1
      value_2=power2

    else if (average_type.eq.'d') then ! average type dB

      value_1=power_dB1
      value_2=power_dB2

    end if
	 
    error=abs(value_1-value_2)
    
    local_ielf_value=local_ielf_value+	&
      error*( (log(function_of_frequency(1)%frequency(sample+1))-log(function_of_frequency(1)%frequency(sample-1)) )/2d0) 	&
            /( log(function_of_frequency(1)%frequency(sample_max+1))-log(function_of_frequency(1)%frequency(sample_min-1)) )

  end do ! next sample
       
  if (local_ielf_value.lt.1d0) then
    judgement='Excellent'
  else if ( (local_ielf_value.ge.1d0).AND.(local_ielf_value.lt.2d0) ) then
    judgement='Very Good'
  else if ( (local_ielf_value.ge.2d0).AND.(local_ielf_value.lt.4d0) ) then
    judgement='Good'
  else if ( (local_ielf_value.ge.4d0).AND.(local_ielf_value.lt.7d0) ) then
    judgement='Moderate'
  else if ( (local_ielf_value.ge.7d0).AND.(local_ielf_value.lt.10d0) ) then
    judgement='Poor'
  else if (local_ielf_value.ge.10d0) then
    judgement='Bad'
  end if
       
  write(*,*)'IELF value:',local_ielf_value,judgement

  write(*,*)'Enter the filename for the IELF data'
  read(*,*)filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
  
  write(local_file_unit,*)'IELF value:',local_ielf_value,judgement
  
  CLOSE(unit=local_file_unit)

  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE IELF
