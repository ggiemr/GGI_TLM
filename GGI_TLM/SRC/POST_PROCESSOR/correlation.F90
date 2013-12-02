! SUBROUTINE correlation_time
! SUBROUTINE correlation_frequency
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
!    correlation_time
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 26/09/2013 CJS
!
!
SUBROUTINE correlation_time

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  integer,parameter	:: max_files=100

  character(len=256)	:: filename

  integer		:: function_number
  
  real*8,allocatable	:: corr(:,:)
  real*8,allocatable	:: mu(:)
  
  integer 		:: i,j
  
  integer 		:: n_samples,sample
  real*8		:: mu1,mu2,ss1,ss2,R12

! START
  
  write(*,*)'Enter the number of time domain quantities for correlation calculation'
  read(*,*)n_functions_of_time
  write(record_user_inputs_unit,*)n_functions_of_time,' number of time domain quantities for correlation calculation '

  if (n_functions_of_time.gt.max_files) then
    write(*,*)'Maximum number of time domain functions exceeded'
    write(*,*)'max_files=',max_files
    write(*,*)'number of time domain quantities for correlation calculation=',n_functions_of_time
    STOP
  end if

  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  do function_number=1,n_functions_of_time
  
    CALL read_Time_Domain_Data(function_number) ! read function of time from output dataset
  
  end do
  
  ALLOCATE( corr(1:n_functions_of_time,1:n_functions_of_time) )  
  ALLOCATE( mu(1:n_functions_of_time) )  

! Check that we have the same number of samples in all functions  

  n_samples=function_of_time(1)%n_timesteps
  
  do i=2,n_functions_of_time
    if ( function_of_time(i)%n_timesteps.ne.n_samples ) then
    
      write(*,*)'Not all the functions have the same number of samples'
      write(*,*)'n_samples (file 1)=',n_samples
      write(*,*)'File number:',i
      write(*,*)'n_samples =',function_of_time(i)%n_timesteps
      STOP
    
    end if
  end do
  
! Calculate mean of each of the functions

  do i=1,n_functions_of_time
    mu(i)=0d0	 
    do sample=1,n_samples
      mu(i)=mu(i)+function_of_time(i)%value(sample)
    end do
    mu(i)=mu(i)/n_samples
  end do
  
! Calculate the correlation matrix

  do i=1,n_functions_of_time
  
    corr(i,i)=1d0   ! autocorrelation =1
    
    do j=i,n_functions_of_time
       
      R12=0D0       
      ss1=0D0       
      ss2=0D0       
  
      do sample=1,n_samples
        R12=R12+((function_of_time(i)%value(sample)-mu(i))*(function_of_time(j)%value(sample)-mu(j)))
        ss1=ss1+(function_of_time(i)%value(sample)-mu(i))**2
        ss2=ss2+(function_of_time(j)%value(sample)-mu(j))**2
      end do

      R12=R12/sqrt(ss1*ss2)
      corr(i,j)=R12
      corr(j,i)=R12
      
    end do
    
  end do
  
  write(*,*)'Enter the filename for the correlation matrix'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)

  OPEN(unit=local_file_unit,file=filename)
  	 
  do i=1,n_functions_of_time
  
    write(local_file_unit,8000)(corr(i,j),j=1,n_functions_of_time)
8000 format(100F10.4)

  end do
	 
  CLOSE(unit=local_file_unit)
  
  DEALLOCATE( corr )  
  DEALLOCATE( mu )  
  
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE correlation_time
!
! NAME
!    correlation_frequency
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 26/09/2013 CJS
!
!
SUBROUTINE correlation_frequency

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  integer,parameter	:: max_files=100

  character(len=256)	:: filename
  character		:: ch
  integer		:: n_functions_of_frequency_in 
  integer		:: function_number
  
  complex*16,allocatable	:: corr(:,:)
  
  integer 		:: i,j
  integer		:: type
  complex*16		:: value,value1,value2
  
  integer 		:: n_samples,sample
  complex*16		:: ss1,ss2,R12
  
  integer		:: n_average

! START
  
  write(*,*)'Enter the number of frequency domain quantities for correlation calculation'
  read(*,*)n_functions_of_frequency_in
  write(record_user_inputs_unit,*)n_functions_of_frequency_in,' number of frequency domain quantities for correlation calculation'

  if (n_functions_of_time.gt.max_files) then
    write(*,*)'Maximum number of frequency domain functions exceeded'
    write(*,*)'max_files=',max_files
    write(*,*)'number of frequency domain quantities for correlation calculation=',n_functions_of_frequency
    STOP
  end if
  
  n_functions_of_frequency=n_functions_of_frequency_in
  n_functions_of_time=0
  
  CALL Allocate_post_data()
  
  do function_number=1,n_functions_of_frequency_in
  
    CALL read_Frequency_Domain_Data(function_number) ! read function of frequency from output dataset
  
  end do
  
  ALLOCATE( corr(1:n_functions_of_frequency_in,1:n_functions_of_frequency_in) )  

! Check that we have the same number of samples in all functions  

  n_samples=function_of_frequency(1)%n_frequencies
  
  do i=2,n_functions_of_frequency_in
    if ( function_of_frequency(i)%n_frequencies.ne.n_samples ) then
    
      write(*,*)'Not all the functions have the same number of samples'
      write(*,*)'n_samples (file 1)=',n_samples
      write(*,*)'File number:',i
      write(*,*)'n_samples =',function_of_frequency(i)%n_frequencies
      STOP
    
    end if
  end do
  
  write(*,*)'Enter the quantity to operate on: Complex, Real, Imaginary or Magnitude'
  read(*,'(A)')ch
  
  if ( (ch.eq.'c').OR.(ch.eq.'C') ) then
  
    write(record_user_inputs_unit,*)'Complex'
    type=0
  
  else if ( (ch.eq.'r').OR.(ch.eq.'R') ) then
  
    write(record_user_inputs_unit,*)'Real'
    type=1
  
  else if ( (ch.eq.'i').OR.(ch.eq.'I') ) then
  
    write(record_user_inputs_unit,*)'Imaginary'
    type=2
  
  else if ( (ch.eq.'m').OR.(ch.eq.'M') ) then
  
    write(record_user_inputs_unit,*)'Magnitude'
    type=3
    
  else
  
    write(*,*)'Quantity to operate on should be one of: Complex, Real, Imaginary or Magnitude'
    STOP
    
  end if
  
! Calculate the correlation matrix

  do i=1,n_functions_of_frequency
    
    do j=i,n_functions_of_frequency
       
      R12=(0D0,0D0)    
      ss1=(0D0,0D0)      
      ss2=(0D0,0D0)      
  
      do sample=1,n_samples
    
        if (type.eq.0) then
          value1=function_of_frequency(i)%value(sample)
          value2=function_of_frequency(j)%value(sample)
        else if (type.eq.1) then
          value1=Dble(function_of_frequency(i)%value(sample))
          value2=Dble(function_of_frequency(j)%value(sample))
        else if (type.eq.2) then
          value1=imag(function_of_frequency(i)%value(sample))
          value2=imag(function_of_frequency(j)%value(sample))
        else if (type.eq.3) then
          value1=function_of_frequency(i)%magnitude(sample)
          value2=function_of_frequency(j)%magnitude(sample)
        end if
      
        R12=R12+value1*conjg(value2) 
        ss1=ss1+value1*conjg(value1)
        ss2=ss2+value2*conjg(value2)
	
      end do

      R12=R12/sqrt(ss1*ss2)
      corr(i,j)=R12
      corr(j,i)=R12
      
    end do
    
  end do
  
  write(*,*)'Enter the filename for the correlation matrix'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)

  OPEN(unit=local_file_unit,file=filename)

  write(local_file_unit,*)'Real part of correlation matrix'
  do i=1,n_functions_of_frequency  
    write(local_file_unit,8000)(Real(corr(i,j)),j=1,n_functions_of_frequency)
8000 format(100F10.4)
  end do

  write(local_file_unit,*)''
  write(local_file_unit,*)'Imaginary part of correlation matrix'
  do i=1,n_functions_of_frequency  
    write(local_file_unit,8000)(Imag(corr(i,j)),j=1,n_functions_of_frequency)
  end do
	 
  CLOSE(unit=local_file_unit)
  
  DEALLOCATE( corr )  
  
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE correlation_frequency
