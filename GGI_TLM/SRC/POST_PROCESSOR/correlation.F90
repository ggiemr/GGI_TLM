! SUBROUTINE correlation_time
! SUBROUTINE correlation_frequency
! SUBROUTINE correlation_function_time
! SUBROUTINE correlation_function_frequency
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
!
! NAME
!    correlation_function_time
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/01/2014 CJS
!
!
SUBROUTINE correlation_function_time

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: filename

  integer		:: function_number
    
  integer 		:: i,j
  
  integer 		:: n_samples1,n_samples2,n_samples3,sample
  real*8		:: dt1,dt2,dt,t1,t2
  
  integer 		:: n_timesteps
  
  real*8		:: min_lag,max_lag
  integer		:: sample_min,sample_max,sample_count
  
  real*8		:: min_tint,max_tint
  integer		:: sample_tint_min,sample_tint_max,n_tint
  
  character 		:: ch
  logical 		:: periodic_extension
  logical 		:: normalise
  
  integer		:: sample1,sample2
  real*8		:: mu1,mu2,f1,f2,CYY,CXX,CXY

! START
  
  n_functions_of_time=3 ! two input functions and one output function

  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  do function_number=1,2
  
    CALL read_Time_Domain_Data(function_number) ! read functions of time 
  
  end do
  
  n_samples1=function_of_time(1)%n_timesteps
  n_samples2=function_of_time(2)%n_timesteps
  dt1=function_of_time(1)%time(2)-function_of_time(1)%time(1)
  dt2=function_of_time(2)%time(2)-function_of_time(2)%time(1)
  
  write(*,*)'Number of time domain samples in file 1=',n_samples1
  write(*,*)'Timestep for data in  file 1           =',dt1
  write(*,*)'Minimum time for data in  file 1       =',function_of_time(1)%time(1)
  write(*,*)'Maximum time for data in  file 1       =',function_of_time(1)%time(n_samples1)
  write(*,*)'Number of time domain samples in file 2=',n_samples2
  write(*,*)'Timestep for data in  file 2           =',dt2
  write(*,*)'Minimum time for data in  file 2       =',function_of_time(2)%time(1)
  write(*,*)'Maximum time for data in  file 2       =',function_of_time(2)%time(n_samples1)

  if (n_samples1.ne.n_samples2) then
    write(*,*)'Error in correlation_function_time'
    write(*,*)'The number of timesteps in the two files should be the same at the moment...'
    STOP
  end if
  
  if (dt1.ne.dt2) then
    write(*,*)'Error in correlation_function_time'
    write(*,*)'The timestep in the two files should be the same...'
    STOP
  end if
  
  if (function_of_time(1)%time(1).ne.function_of_time(2)%time(1)) then
    write(*,*)'Error in correlation_function_time'
    write(*,*)'The time origin in the two files should be the same...'
    STOP
  end if
  
  dt=dt1
  
  write(*,*)'Do you want to assume a periodic extension of the two input functions? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  periodic_extension=.FALSE.
  if (ch.eq.'y') then
    periodic_extension=.TRUE.
  end if

  write(*,*)'Enter the minimum time lag for the correlation function, (zero for min value)'
  read(*,*)min_lag
  write(record_user_inputs_unit,*)min_lag,' :minimum time lag for the correlation function, (zero for min value)'

  write(*,*)'Enter the maximum time lag for the correlation function, (zero for max value)'
  read(*,*)max_lag
  write(record_user_inputs_unit,*)max_lag,' :maximum time lag for the correlation function, (zero for max value)'
  
  if (min_lag.eq.0d0) then
    if (.NOT.periodic_extension) then
      min_lag=-function_of_time(1)%time(n_samples1)
    else
      min_lag=0d0
    end if
  end if
  
  if (max_lag.eq.0d0) then
    if (.NOT.periodic_extension) then
      max_lag=function_of_time(1)%time(n_samples1)
    else
      max_lag=function_of_time(1)%time(n_samples1)
    end if
  end if
  
  sample_min=NINT(min_lag/dt)
  sample_max=NINT(max_lag/dt)
  write(*,*)'Minimum lag time=',min_lag,' sample number',sample_min
  write(*,*)'Maximum lag time=',max_lag,' sample number',sample_max

  n_timesteps=(sample_max-sample_min)+1
  write(*,*)'Number of time domain samples in time domain correlation file=',n_timesteps

  write(*,*)'Enter the minimum integration time for the correlation function, (zero for min value)'
  read(*,*)min_tint
  write(record_user_inputs_unit,*)min_tint,' :minimum integration time for the correlation function, (zero for min value)'

  write(*,*)'Enter the maximum integration time for the correlation function, (zero for max value)'
  read(*,*)max_tint
  write(record_user_inputs_unit,*)max_tint,' :maximum integration time for the correlation function, (zero for max value)'
  
  if (min_tint.eq.0d0) then
    min_tint=function_of_time(1)%time(1)
  end if
  
  if (max_tint.eq.0d0) then
    max_tint=function_of_time(1)%time(n_samples1)
  end if

! work out the integration range  
  sample_tint_min=1
  sample_tint_max=n_samples1
  
  do sample=1,n_samples1-1
  
    t1=function_of_time(1)%time(sample)
    t2=function_of_time(1)%time(sample+1)
    
    if ( (t1.LE.min_tint).AND.(t2.GT.min_tint) ) sample_tint_min=sample
    if ( (t1.LT.max_tint).AND.(t2.GE.max_tint) ) sample_tint_max=sample+1
    
  end do
  
  n_tint=sample_tint_max-sample_tint_min+1
  
  write(*,*)'Minimum integration time=',min_tint,' sample number',sample_tint_min
  write(*,*)'Maximum integration time=',max_tint,' sample number',sample_tint_max
  write(*,*)'Number of integration samples=',n_tint
  
  write(*,*)'Do you want to normalise the correlation function ( C=CXY/sqrt(CXX*CYY) )? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  normalise=.FALSE.
  if (ch.eq.'y') then
    normalise=.TRUE.
  end if

! allocate time domain output function  
  
  function_of_time(3)%n_timesteps=n_timesteps
  ALLOCATE ( function_of_time(3)%time(1:n_timesteps) )
  ALLOCATE ( function_of_time(3)%value(1:n_timesteps) )
  
! Calculate mean and variance of each of the functions over the integration range

  mu1=0d0      
  do sample=sample_tint_min,sample_tint_max
    mu1=mu1+function_of_time(1)%value(sample)
  end do
  mu1=mu1/n_tint
  CXX=0d0
  do sample=sample_tint_min,sample_tint_max
    CXX=CXX+(function_of_time(1)%value(sample)-mu1)*(function_of_time(1)%value(sample)-mu1)*dt
  end do
!  CXX=CXX/(n_tint*dt)
  
  write(*,*)'Function 1: mean,',mu1,' variance',CXX
  
  mu2=0d0      
  do sample=sample_tint_min,sample_tint_max
    mu2=mu2+function_of_time(2)%value(sample)
  end do
  mu2=mu2/n_tint
  CYY=0d0
  do sample=sample_tint_min,sample_tint_max
    CYY=CYY+(function_of_time(2)%value(sample)-mu2)*(function_of_time(2)%value(sample)-mu2)*dt
  end do
!  CYY=CYY/(n_tint*dt)
  
  write(*,*)'Function 2: mean,',mu2,' variance',CYY
  
! Calculate the correlation function

  sample_count=0

  do i=sample_min,sample_max  ! time lag loop
  
    sample_count=sample_count+1
    function_of_time(3)%time(sample_count)=(sample_min+sample_count-1)*dt
    
    CXY=0d0
    
    if (periodic_extension) then
    
      do j=1,n_tint  ! time integration loop
          
        f1=0d0
        f2=0d0
	
        sample1=mod(j-1,n_tint)
        if (sample1.lt.0) sample1=n_tint+sample1
        sample1=sample1+sample_tint_min

        sample2=mod(j-i-1,n_tint)
        if (sample2.lt.0) sample2=n_tint+sample2
        sample2=sample2+sample_tint_min
      
        if ((sample1.GE.sample_tint_min).AND.(sample1.LE.sample_tint_max)) then
	  f1=function_of_time(1)%value(sample1)-mu1
	else
	  write(*,*)'Error: sample1=',sample1,' range:',sample_tint_min,sample_tint_max
	end if
        if ((sample2.GE.sample_tint_min).AND.(sample2.LE.sample_tint_max)) then
	  f2=function_of_time(2)%value(sample2)-mu2
	else
	  write(*,*)'Error: sample2=',sample2,' range:',sample_tint_min,sample_tint_max
	end if
      
        CXY=CXY+f1*f2*dt
      
      end do  ! next point in time integration
      
!      CXY=CXY/(n_tint*dt)
 
    else
! no periodic extension, assume the value is zero outside the specified range 
    
      do j=sample_tint_min,sample_tint_max  ! time integration loop
          
        f1=0d0
        f2=0d0
        sample1=j
        sample2=j-i
        if ((sample1.GE.sample_tint_min).AND.(sample1.LE.sample_tint_max)) f1=function_of_time(1)%value(sample1)-mu1
        if ((sample2.GE.sample_tint_min).AND.(sample2.LE.sample_tint_max)) f2=function_of_time(2)%value(sample2)-mu2
      
        CXY=CXY+f1*f2*dt
       
      end do  ! next point in time integration
      
!      CXY=CXY/(n_tint*dt)
         
    end if  ! periodic extension
    
    if (normalise) then
    
      if ( (CXX.NE.0d0).AND.(CYY.NE.0d0) ) then
        function_of_time(3)%value(sample_count)=CXY/sqrt(CXX*CYY)
      else
        function_of_time(3)%value(sample_count)=0d0
      end if
    else
    
      function_of_time(3)%value(sample_count)=CXY
    
    end if ! noralise
    
  end do  ! next correlation time lag
  
! Write the correlation data set to file
  
  CALL write_time_domain_data(3)
  
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE correlation_function_time
!
! NAME
!    correlation_function_frequency
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/01/2014 CJS
!
!
SUBROUTINE correlation_function_frequency

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: filename

  integer		:: function_number
    
  integer 		:: i,j
  
  integer 		:: n_samples1,n_samples2,n_samples3,sample
  real*8		:: df1,df2,df,f1,f2
  
  integer 		:: n_frequencies
  
  real*8		:: min_lag,max_lag
  integer		:: sample_min,sample_max,sample_count
  
  real*8		:: min_fint,max_fint
  integer		:: sample_fint_min,sample_fint_max,n_fint
  
  character 		:: ch
  logical 		:: periodic_extension
  logical 		:: normalise
  
  integer		:: sample1,sample2
  complex*16		:: mu1,mu2,value1,value2,CYY,CXX,CXY

! START
  
  n_functions_of_time=0 ! two input functions and one output function

  n_functions_of_frequency=3 ! two input functions and one output function
  
  CALL Allocate_post_data()
  
  do function_number=1,2
  
    CALL read_Frequency_Domain_Data(function_number) ! read functions of frequency
  
  end do
  
  n_samples1=function_of_frequency(1)%n_frequencies
  n_samples2=function_of_frequency(2)%n_frequencies
  df1=function_of_frequency(1)%frequency(2)-function_of_frequency(1)%frequency(1)
  df2=function_of_frequency(2)%frequency(2)-function_of_frequency(2)%frequency(1)
  
  write(*,*)'Number of frequency domain samples in file 1=',n_samples1
  write(*,*)'Frequency step for data in  file 1           =',df1
  write(*,*)'Minimum frequency for data in  file 1       =',function_of_frequency(1)%frequency(1)
  write(*,*)'Maximum frequency for data in  file 1       =',function_of_frequency(1)%frequency(n_samples1)
  write(*,*)'Number of frequency domain samples in file 2=',n_samples2
  write(*,*)'Frequency step for data in  file 2           =',df2
  write(*,*)'Minimum frequency for data in  file 2       =',function_of_frequency(2)%frequency(1)
  write(*,*)'Maximum frequency for data in  file 2       =',function_of_frequency(2)%frequency(n_samples1)

  if (n_samples1.ne.n_samples2) then
    write(*,*)'Error in correlation_function_frequency'
    write(*,*)'The number of frequency steps in the two files should be the same at the moment...'
    STOP
  end if
  
  if (df1.ne.df2) then
    write(*,*)'Error in correlation_function_frequency'
    write(*,*)'The frequency step in the two files should be the same...'
    STOP
  end if
  
  if (function_of_frequency(1)%frequency(1).ne.function_of_frequency(2)%frequency(1)) then
    write(*,*)'Error in correlation_function_frequency'
    write(*,*)'The frequency origin in the two files should be the same...'
    STOP
  end if
  
  df=df1
  
  write(*,*)'Do you want to assume a periodic extension of the two input functions? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  periodic_extension=.FALSE.
  if (ch.eq.'y') then
    periodic_extension=.TRUE.
  end if

  write(*,*)'Enter the minimum frequency lag for the correlation function, (zero for min value)'
  read(*,*)min_lag
  write(record_user_inputs_unit,*)min_lag,' :minimum frequency lag for the correlation function, (zero for min value)'

  write(*,*)'Enter the maximum frequency lag for the correlation function, (zero for max value)'
  read(*,*)max_lag
  write(record_user_inputs_unit,*)max_lag,' :maximum frequency lag for the correlation function, (zero for max value)'

  if (min_lag.eq.0d0) then
    min_lag=function_of_frequency(1)%frequency(1)
  end if
  
  if (max_lag.eq.0d0) then
    max_lag=function_of_frequency(1)%frequency(n_samples1)
  end if
  
  sample_min=NINT(min_lag/df)
  sample_max=NINT(max_lag/df)
  write(*,*)'Minimum lag frequency=',min_lag,' sample number',sample_min
  write(*,*)'Maximum lag frequency=',max_lag,' sample number',sample_max

  n_frequencies=(sample_max-sample_min)+1
  write(*,*)'Number of frequency domain samples in frequency domain correlation file=',n_frequencies

  write(*,*)'Enter the minimum integration frequency for the correlation function, (zero for min value)'
  read(*,*)min_fint
  write(record_user_inputs_unit,*)min_fint,' :minimum integration frequency for the correlation function, (zero for min value)'

  write(*,*)'Enter the maximum integration frequency for the correlation function, (zero for max value)'
  read(*,*)max_fint
  write(record_user_inputs_unit,*)max_fint,' :maximum integration frequency for the correlation function, (zero for max value)'
  
  if (min_fint.eq.0d0) then
    min_fint=function_of_frequency(1)%frequency(1)
  end if
  
  if (max_fint.eq.0d0) then
    max_fint=function_of_frequency(1)%frequency(n_samples1)
  end if

! work out the integration range  
  sample_fint_min=1
  sample_fint_max=n_samples1
  
  do sample=1,n_samples1-1
  
    f1=function_of_frequency(1)%frequency(sample)
    f2=function_of_frequency(1)%frequency(sample+1)
    
    if ( (f1.LE.min_fint).AND.(f2.GT.min_fint) ) sample_fint_min=sample
    if ( (f1.LT.max_fint).AND.(f2.GE.max_fint) ) sample_fint_max=sample+1
    
  end do
  
  n_fint=sample_fint_max-sample_fint_min+1
  
  write(*,*)'Minimum integration frequency=',min_fint,' sample number',sample_fint_min
  write(*,*)'Maximum integration frequency=',max_fint,' sample number',sample_fint_max
  write(*,*)'Number of integration samples=',n_fint
  
  write(*,*)'Do you want to normalise the correlation function ( C=CXY/sqrt(CXX*CYY) )? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  normalise=.FALSE.
  if (ch.eq.'y') then
    normalise=.TRUE.
  end if

! allocate frequency domain output function  
  
  function_of_frequency(3)%n_frequencies=n_frequencies
  ALLOCATE ( function_of_frequency(3)%frequency(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(3)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(3)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(3)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(3)%dB(1:n_frequencies) )
  
  function_of_frequency(3)%frequency(1:n_frequencies)=0d0
  function_of_frequency(3)%value(1:n_frequencies)=(0d0,0d0)
  function_of_frequency(3)%magnitude(1:n_frequencies)=0d0
  function_of_frequency(3)%phase(1:n_frequencies)=0d0
  function_of_frequency(3)%dB(1:n_frequencies)=0d0
  
! Calculate mean and variance of each of the functions over the integration range

  mu1=0d0      
  do sample=sample_fint_min,sample_fint_max
    mu1=mu1+function_of_frequency(1)%value(sample)
  end do
  mu1=mu1/n_fint
  CXX=0d0
  do sample=sample_fint_min,sample_fint_max
    CXX=CXX+(function_of_frequency(1)%value(sample)-mu1)*(function_of_frequency(1)%value(sample)-mu1)
  end do
  CXX=CXX/n_fint
  
  write(*,*)'Function 1: mean    :',mu1
  write(*,*)'Function 1: variance:',CXX
  
  mu2=0d0      
  do sample=sample_fint_min,sample_fint_max
    mu2=mu2+function_of_frequency(2)%value(sample)
  end do
  mu2=mu2/n_fint
  CYY=0d0
  do sample=sample_fint_min,sample_fint_max
    CYY=CYY+(function_of_frequency(2)%value(sample)-mu2)*(function_of_frequency(2)%value(sample)-mu2)
  end do
  CYY=CYY/n_fint
  
  write(*,*)'Function 2: mean    :',mu2
  write(*,*)'Function 2: variance:',CYY
  
! Calculate the correlation function

  sample_count=0

  do i=sample_min,sample_max  ! frequency lag loop
  
    sample_count=sample_count+1
    function_of_frequency(3)%frequency(sample_count)=(sample_min+sample_count-1)*df
    
    CXY=0d0
    
    if (periodic_extension) then
    
      do j=1,n_fint  ! frequency integration loop
          
        value1=0d0
        value2=0d0
	
        sample1=mod(j-1,n_fint)
        if (sample1.lt.0) sample1=n_fint+sample1
        sample1=sample1+sample_fint_min

        sample2=mod(j-i-1,n_fint)
        if (sample2.lt.0) sample2=n_fint+sample2
        sample2=sample2+sample_fint_min
      
        if ((sample1.GE.sample_fint_min).AND.(sample1.LE.sample_fint_max)) then
	  value1=function_of_frequency(1)%value(sample1)-mu1
	else
	  write(*,*)'Error: sample1=',sample1,' range:',sample_fint_min,sample_fint_max
	end if
        if ((sample2.GE.sample_fint_min).AND.(sample2.LE.sample_fint_max)) then
	  value2=function_of_frequency(2)%value(sample2)-mu2
	else
	  write(*,*)'Error: sample2=',sample2,' range:',sample_fint_min,sample_fint_max
	end if
      
        CXY=CXY+value1*value2
      
      end do  ! next point in frequency integration
      
      CXY=CXY/n_fint
 
    else
! no periodic extension, assume the value is zero outside the specified range 
    
      do j=sample_fint_min,sample_fint_max  ! frequency integration loop
          
        value1=0d0
        value2=0d0
        sample1=j
        sample2=j-i
        if ((sample1.GE.sample_fint_min).AND.(sample1.LE.sample_fint_max)) value1=function_of_frequency(1)%value(sample1)-mu1
        if ((sample2.GE.sample_fint_min).AND.(sample2.LE.sample_fint_max)) value2=function_of_frequency(2)%value(sample2)-mu2
      
        CXY=CXY+value1*value2
       
      end do  ! next point in frequency integration
      
      CXY=CXY/n_fint
         
    end if  ! periodic extension
    
    if (normalise) then
    
      if ( (CXX.NE.0d0).AND.(CYY.NE.0d0) ) then
        function_of_frequency(3)%value(sample_count)=CXY/sqrt(CXX*CYY)
      else
        function_of_frequency(3)%value(sample_count)=0d0
      end if
    else
    
      function_of_frequency(3)%value(sample_count)=CXY
    
    end if ! noralise
    
    function_of_frequency(3)%magnitude(sample_count)=	&
                    abs(function_of_frequency(3)%value(sample_count))
    function_of_frequency(3)%phase(sample_count)=	&
                    atan2( imag(function_of_frequency(3)%value(sample_count)), &
                           dble(function_of_frequency(3)%value(sample_count))   )
    function_of_frequency(3)%dB(sample_count)=	&
                    20d0*log10(function_of_frequency(3)%magnitude(sample_count))
    
  end do  ! next correlation frequency lag
  
! Write the correlation data set to file
  
  CALL write_frequency_domain_data(3)
  
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE correlation_function_frequency
