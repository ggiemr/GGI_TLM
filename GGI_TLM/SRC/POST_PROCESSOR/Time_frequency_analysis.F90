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
! SUBROUTINE Time_frequency_analysis
!
! NAME
!    Time_frequency_analysis
!
! DESCRIPTION
!    Time-frequency analysis i.e. perform a short time FFT on a windowed version of the input time sequence,
!    allowing the window to move through the input time domain data
!     
! COMMENTS
!     Note that this process operates on simple data files with no header data i.e. it could come from measurement or wherever...
!     The input format is assumed to be two column data:
!
!     time    value
!
!     GGI_TLM output files will require format conversion before using this process...
!
! HISTORY
!
!     started 19/11/2013 CJS
!     10/3/2016 CJS. Try to make the system work better for very large datasets
!
SUBROUTINE Time_frequency_analysis

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  character(len=256)	:: filename
 
  logical		:: file_exists
  
  real*8	:: t_in,value_in
  
  real*8	:: fmin,fmax,fstep
  integer	:: n_frequencies_in_range,n_frequencies,equivalent_n_frequencies_in_full_range
  integer	:: frequency_loop
  integer	:: timestep,n_timesteps,n_timesteps_read
  
  real*8	:: dt
  
  real*8,allocatable		:: time(:)
  real*8,allocatable		:: value(:)
  complex*16,allocatable	:: x(:)
  real*8,allocatable		:: frequency(:)
  
  integer	:: window_width,WW2,window_width_applied
  integer 	:: n_time_windows,time_window
  integer	:: window_centre_min,window_centre_max
  integer	:: window_centre,window_min,window_max
  real*8	:: window_centre_step_real
  real*8	:: window_centre_real
  real*8	:: time_window_centre
  
  real*8	:: fmin_out,fmax_out
  
  integer	:: min_nfreq

  character(len=256)	:: command
  
  integer 	:: i,loop,count
  real*8	:: fstep_out
  integer	:: ifstride
  
  real*8        :: sumsqr
  character     :: ch
  logical	:: plot_dB
  
! START

! read the input data set
  write(*,*)'Output files:'
  
  command='ls -ltr '
  CALL system(command)

  write(*,*)'Enter the time domain filename'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    STOP
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  write(post_process_info_unit,*)'	Time domain data filename:',trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
    
  write(*,*)'Opened file:',trim(filename)
       
  do loop=1,2
  
  timestep=0
    
10  CONTINUE

      read(local_file_unit,*,end=1000)t_in,value_in
            
      timestep=timestep+1
	
      if (loop.eq.2) then

        time(timestep)=t_in
        value(timestep)=value_in
        
      end if ! loop.EQ.2
	    
    GOTO 10  ! read next line of data
   
1000 CONTINUE

    n_timesteps_read=timestep
    
    if (loop.eq.1) then
      ALLOCATE ( time(1:n_timesteps_read) )
      ALLOCATE ( value(1:n_timesteps_read) )
    end if
    
    rewind(unit=local_file_unit)
    
  end do ! next loop
  	 
  CLOSE(unit=local_file_unit)

  n_timesteps=n_timesteps_read
  dt=time(2)-time(1)
  
! Write a summary of the time domain data

  write(*,*)'__________________________________________________'
  write(*,*)''
  write(*,*)'First time value:',time(1)
  write(*,*)'Last time value :',time(n_timesteps)
  write(*,*)'Number of timesteps :',n_timesteps
  write(*,*)'timestep :',dt
  write(*,*)''
  
! read the time window width and the number of time windows to analyse

  write(*,*)'Enter the width of the time window (in timesteps) for time-frequency analysis'
  read(*,*)window_width
  write(record_user_inputs_unit,*)window_width,' Window width'
 
  write(*,*)'Enter the number of time windows to analyse for time-frequency analysis'
  read(*,*)n_time_windows
  write(record_user_inputs_unit,*)n_time_windows,' Number of time windows'

! work out the frequency range of the FFT output
  fmin=0d0
  fmax=1/dt
  
  write(*,*)'FFT Fmin=',fmin,' Hz'
  write(*,*)'FFT Fmax=',fmax,' Hz'
 
  write(*,*)'Enter the minimum frequency to ouptut:'
  read(*,*)fmin_out
  write(record_user_inputs_unit,*)fmin_out,'   Minimum frequency to ouptut'
 
  write(*,*)'Enter the maximum frequency to ouptut:'
  read(*,*)fmax_out
  write(record_user_inputs_unit,*)fmax_out,'   Maximum frequency to ouptut'
 
  write(*,*)'Enter the approximate number of frequency samples in this range to ouptut:'
  read(*,*)n_frequencies_in_range
  write(record_user_inputs_unit,*)n_frequencies_in_range,'   Approximate number of frequency samples in this range to ouptut'

! work out the equivalent number of frequency samples in the full FFT frequency range
! if we had asked for this frequency resolution
  equivalent_n_frequencies_in_full_range=NINT(dble(n_frequencies_in_range)*fmax/(fmax_out-fmin_out))
  
  write(*,*)'equivalent number of frequency samples in the full FFT frequency range:',equivalent_n_frequencies_in_full_range
  
! set the minimum number of frequencies to be equal to the next power of 2 above the equivalent number of frequency samples in the full FFT
  min_nfreq=max(equivalent_n_frequencies_in_full_range,window_width)
  
  i=1
20  CONTINUE
    i=i*2
    if (i.ge.min_nfreq) then
      n_frequencies=i
    else
      GOTO 20
    end if

  write(*,*)'Number of samples required in FFT=',n_frequencies

  fmin=0d0
  fstep=1d0/(n_frequencies*dt)
  fmax=fstep*(n_frequencies-1)

  write(*,*)'FFT Fmin=',fmin,' Hz'
  write(*,*)'FFT Fmax=',fmax,' Hz'
  write(*,*)'FFT delta f=',fstep,' Hz'
  write(*,*)'Number of samples in FFT=',n_frequencies

! This ensures that we don't output too many samples if the frequency resolution from FFT is much smaller than
! that requested  
  fstep_out=(fmax_out-fmin_out)/n_frequencies_in_range
  ifstride=NINT(fstep_out/fstep)
  if (ifstride.LT.1) ifstride=1

  write(*,*)'frequency output step (equivalent_n_frequencies_in_full_range/n_frequencies_in_range)',ifstride
  
  plot_dB=.FALSE.
  write(*,*)'plot in dB? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch
  CALL convert_to_lower_case(ch,1)
  if (ch.eq.'y') then
    plot_dB=.TRUE.
  end if

! Open the output file 
  write(*,*)'Enter the time frequency output filename'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)

! setup the time window loop

  WW2=window_width/2
  window_centre_min=WW2
  window_centre_max=n_timesteps-WW2
  
  if (n_time_windows.GT.1) then
    window_centre_step_real=dble(window_centre_max-window_centre_min)/dble(n_time_windows-1)
  else
    window_centre_step_real=0d0
  end if
  
  write(*,*)'Window width=',window_width

! Allocate memory for FFT and Generate the FFT frequency list
  
  ALLOCATE( frequency(1:n_frequencies) )
  ALLOCATE( x(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
      
    frequency(frequency_loop)=fmin+(frequency_loop-1)*fstep
  
  end do
  
! loop over the time window

  do time_window=1,n_time_windows
  
    window_centre_real=dble(window_centre_min)+dble(time_window-1)*window_centre_step_real
    window_centre=NINT(window_centre_real)
    window_min=window_centre-WW2
    window_max=window_centre+WW2
    window_width_applied=window_max-window_min+1
    time_window_centre=(window_centre_real-1d0)*dt
    
! ensure that the window doesn't go outside the bounds of the time array
    if ( (window_min.LT.1).OR.(window_max.GT.n_timesteps) ) then
      write(*,*)'Error: window out of range:',window_min,window_max,' max=',n_frequencies
    end if
    
    window_min=max(1,window_min)
    window_max=min(n_timesteps,window_max)
    
! copy the windowed time domain data into the temporary complex array, x

    x(:)=(0d0,0d0)
  
    x(1:window_width_applied)=dcmplx(value(window_min:window_max))

! Fourier transform using the FFT algorithm  
  
    CALL FFT(x,n_frequencies)
  
! write the required FFT data into the output file
    do frequency_loop=1,n_frequencies,ifstride
  
      if ( (frequency(frequency_loop).GE.fmin_out).AND.(frequency(frequency_loop).LE.fmax_out) ) then
    
        sumsqr=0d0
        count=0
        do i=frequency_loop,min(frequency_loop+ifstride,n_frequencies)
          count=count+1
          sumsqr=sumsqr+(abs(x(i)))**2
        end do
        if (count.NE.0) then
          sumsqr=sumsqr/count
        end if

! OLD        
!        write(local_file_unit,8000)frequency(frequency_loop),time_window_centre,dble(x(frequency_loop)),	&
!                                   imag(x(frequency_loop)),abs(x(frequency_loop))

        if (plot_dB) then
          write(local_file_unit,8000)frequency(frequency_loop),time_window_centre,10d0*log10(sumsqr)
        else
          write(local_file_unit,8000)frequency(frequency_loop),time_window_centre,sqrt(sumsqr)    
        end if
        
8000  format(3E16.6)
      
      end if
    
    end do ! next frequency value
    
  end do ! next time window
  
  CLOSE(unit=local_file_unit)

  DEALLOCATE ( time )
  DEALLOCATE ( value )
  DEALLOCATE( frequency )
  DEALLOCATE( x )

  RETURN
  
END SUBROUTINE Time_frequency_analysis
