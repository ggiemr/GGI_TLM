SUBROUTINE convert_data_to_sound

! read time domain data, apply some data processing (filters) and write to a .wav file
! i.e. an output format that can be played through a speaker

USE constants
USE file_information

IMPLICIT NONE

integer :: nt,nt_pad
character(LEN=2) :: ch2
character        :: ch
character*256    :: ipfilename
character*256    :: filename
character*256    :: opfilename
real*8           :: tin,ftin
integer :: readloop

real*8 :: dt

integer :: n_head,t_col,V_col,n_col

real*8,allocatable :: data_line(:)

real*8,allocatable :: t(:),ft(:)
real*8 :: tmax

! Functions to be transformed are complex here
real*8,allocatable     :: frequency(:)
complex*16,allocatable :: ft_filtered(:)
complex*16,allocatable :: ft_impulse(:)

integer :: n_filter,filter
integer :: order
real*8  :: fc,wc

integer :: n_sub_segments
integer :: nt_sub,nt_sub2
integer :: step_sub
integer :: n_sub_plot
integer :: nt_sub_pad,nt_sub_pad2
real*8    :: sub_sample_period
real*8    :: dt_sub,t_sub_period
real*8    :: tstart,tstop,time
integer   :: ntstart,istart
real*8    :: window_spacing

real*8  :: t1,t2
integer :: last_sample,i_time
real*8  :: interpolated_value,value_1,value_2

real*8,allocatable :: t_sub(:)
real*8,allocatable :: ft_sub(:)

complex*16 :: Vdc

integer i,ii,isegment
  
logical		:: file_exists

integer :: sample,n_samples

real :: real_max,real_min,real_avg,real_range

character*4 ChunkID
integer*4 ChunkSize
character*4 Format
character*4 Subchunk1ID

integer*4 Subchunk1Size

integer*2 AudioFormat
integer*2 NumChannels
integer*4 SampleRate
integer*4 ByteRate
integer*2 BlockAlign
integer*2 BitsPerSample

character*4 Subchunk2ID
integer*4 Subchunk2Size

integer*2,allocatable ::  int_data(:)

!logical,parameter :: write_temp_files=.FALSE.
logical,parameter :: write_temp_files=.TRUE.

! START

! 1. READ THE INPUT DATASET

write(*,*)'Enter the filename for the time domain data'
read(*,'(A)')ipfilename

inquire(file=trim(ipfilename),exist=file_exists)
if (.NOT.file_exists) then
  write(*,*)'Error file does not exist:',trim(ipfilename)
  STOP
end if
write(record_user_inputs_unit,'(A)')trim(ipfilename)
write(post_process_info_unit,*)'Time domain data filename:',trim(ipfilename)

open(unit=local_file_unit,file=trim(ipfilename),STATUS='OLD')

write(*,*)'Enter the number of header lines in the file'
read(*,*)n_head
write(record_user_inputs_unit,*)n_head,'  # number of header lines in the file'

write(*,*)'Enter the column number for time data'
read(*,*)t_col
write(record_user_inputs_unit,*)t_col,'  # time data column'

write(*,*)'Enter the column number for function data'
read(*,*)V_col
write(record_user_inputs_unit,*)V_col,'  # Voltage data column'

n_col=max(t_col,V_col)

ALLOCATE( data_line(1:n_col) )

do readloop=1,2

! read file header
  do i=1,n_head
    read(local_file_unit,*)
  end do

  nt=0

10 CONTINUE

  read(local_file_unit,*,END=100)(data_line(i),i=1,n_col)
  tin =data_line(t_col)
  ftin=data_line(V_col)
  nt=nt+1
  if (readloop.EQ.2) then
    t(nt)=tin
    ft(nt)=ftin
  end if

  GOTO 10
  
100 CONTINUE

  if (readloop.EQ.1) then

! work out the power of two greater than or equal to nt
    nt_pad=2
    do while (nt_pad.LT.nt)
      nt_pad=nt_pad*2
    end do

! Allocate memory for the input time domain dataset

    ALLOCATE( ft(1:nt_pad) )
    ALLOCATE( t(1:nt_pad) )
    ft(1:nt_pad)=0d0
    t(1:nt_pad) =0d0  
    rewind(local_file_unit)
  end if
  
end do ! next readloop

CLOSE(unit=local_file_unit)
DEALLOCATE( data_line )

dt=(t(nt)-t(1))/dble(nt-1)

write(*,*)
write(*,*)'Number of samples read, nt=',nt
write(*,*)'Time of first sample     =',t(1)
write(*,*)'Time of last sample      =',t(nt)
write(*,*)'Timestep of input signal =',dt

tmax=t(nt)-t(1)

write(*,*)'Time period of the input signal=',tmax
write(*,*)'Rescaling time axis from 0 to tmax, tmax=',tmax
write(*,*)
write(*,*)'Number of samples with zero padding to the next power of 2, nt_pad=',nt_pad
write(*,*)'(Note: zero padding is only used here for the purposes of applying filters)'

! fill the time axis data including the padded values
do i=1,nt_pad
  t(i)=dble(i-1)*dt
end do

write(*,*)'Time of last sample with zero padding=',t(nt_pad)

!2. APPLY LOW PASS AND HIGH PASS FILTERS TO INPUT TIME DOMAIN DATA AS REQUIRED

! Work out the frequency for each sample
ALLOCATE( frequency(1:nt_pad) )
CALL FFT_set_frequencies(frequency,nt_pad,dt)

ALLOCATE( ft_filtered(1:nt_pad) )

write(*,*)'Enter the number of filters to apply to the input time domain data'
read(*,*)n_filter
write(record_user_inputs_unit,*)n_filter,'  # number of filters to apply'
write(*,*)'number of filters=',n_filter

ft_filtered(1:nt_pad)=dcmplx(ft(1:nt_pad))

! allocate data for the impulse response
ALLOCATE( ft_impulse(1:nt_pad) )
ft_impulse(1:nt_pad)=(0d0,0d0)
ft_impulse(1)=(1d0,0d0)

! Fourier Transform the input funtion
  
CALL FFT(ft_impulse,nt_pad)
CALL FFT(ft_filtered,nt_pad)
  
if (write_temp_files) then
  filename='ft_freq.dat'
  CALL write_FFT_data(frequency,ft_filtered,nt_pad,filename,local_file_unit)
end if

do filter=1,n_filter

  write(*,*)'Enter the filter type (LPF or HPF)'
  read(*,*)ch
  write(record_user_inputs_unit,'(A,A)')ch,'  # filter type (LPF or HPF)'
  write(*,*)'Filter type=',ch

  write(*,*)'Enter the filter order (integer)'
  read(*,*)order
  write(record_user_inputs_unit,*)order,'  # filter order'
  write(*,*)'Order=',order
  
  write(*,*)'Enter the filter 3dB frequency (Hz)'
  read(*,*)fc
  write(record_user_inputs_unit,*)fc,'  # filter cutoff frequency'
  write(*,*)'cutoff frequency=',fc
  wc=2d0*pi*fc
  
  if( (ch.eq.'l').OR.(ch.eq.'L') ) then
  
! Apply a low pass filter function in the frequency domain
  
    CALL butterworth_LPF(frequency,ft_impulse,nt_pad,order,fc,dt)
    CALL butterworth_LPF(frequency,ft_filtered,nt_pad,order,fc,dt)
    
  else if( (ch.eq.'h').OR.(ch.eq.'H') ) then
  
! Apply a high pass filter function in the frequency domain
  
    CALL butterworth_HPF(frequency,ft_impulse,nt_pad,order,fc,dt)
    CALL butterworth_HPF(frequency,ft_filtered,nt_pad,order,fc,dt)
   
  else
    write(*,*)'The filter should be LPF or HPF'
    STOP 1
  end if

! Inverse Fourier Transform the filtered data to give the time response
  
end do ! apply next filter function
  
if (write_temp_files) then
  filename='filter_transfer_function.dat'
  CALL write_FFT_data(frequency,ft_impulse,nt_pad,filename,local_file_unit)
  
  filename='ft_filtered.dat'
  CALL write_FFT_data(frequency,ft_filtered,nt_pad,filename,local_file_unit)
end if
 
CALL INVERSE_FFT(ft_impulse,nt_pad)
CALL INVERSE_FFT(ft_filtered,nt_pad)
  
if (write_temp_files) then
  filename='filter_impulse_response.dat'
  CALL write_FFT_data(t,ft_impulse,nt_pad,filename,local_file_unit)
  
  filename='ft_out.dat'
  CALL write_FFT_data(t,ft_filtered,nt_pad,filename,local_file_unit)
end if

!3 EXTRACT SUB-SAMPLE FROM THE INPUT DATASET

write(record_user_inputs_unit,*)'Number of sub-segments of data to convert =1'

write(*,*)'Total number of input samples=',nt

write(*,*)'Period of input signal=',(t(nt)-t(1))

write(*,*)'Enter the start time of the data sub-segment'
read(*,*)tstart
write(record_user_inputs_unit,*)tstart,'  # start time of the data sub-segment'
write(*,*)'start time of the data sub-segment=',tstart

write(*,*)'Enter the end time of the data sub-segment'
read(*,*)tstop
write(record_user_inputs_unit,*)tstop,'  # end time of the data sub-segment'
write(*,*)'end time of the data sub-segment=',tstop

t_sub_period=tstop-tstart

write(*,*)'Enter the number of samples in the data sub-segment'
read(*,*)nt_sub
write(record_user_inputs_unit,*)nt_sub,'  # number of samples in the data sub-segment'
write(*,*)'number of samples in the data sub-segment=',nt_sub

n_sub_plot=1

dt_sub=t_sub_period/dble(nt_sub-1)       ! timestep for the sub-sampled segments

write(*,*)'Timestep for sub-sampled segments =',dt_sub

! set the time axis data for sub-segments incuding zero padding
ALLOCATE( t_sub(1:nt_sub) ) 
do i=1,nt_sub
  t_sub(i)=dble(i-1)*dt_sub
end do

ALLOCATE( ft_sub(1:nt_sub) )
ft_sub(1:nt_sub)=(0d0,0d0)

write(*,*)
write(*,*)'extract the sub-segments from the original filtered dataset'

  tstart=0d0
  
  write(*,*)'Sub segment initial time: ',tstart  
    
! Find the timesteps that bracket the initial time
    
  do ii=1,nt_pad-1
    
    t1=t(ii)
    t2=t(ii+1)
      
    if ( (tstart.GE.t1).AND.(tstart.LE.t2) ) then
! the initial time is bracketed by timesteps ii and ii+1
      GOTO 7000
    end if
    
  end do
    
  write(*,*)'ERROR: Unable to find the intiial timestep for this sub-segment data'
  STOP 1
    
7000 CONTINUE   ! Jump here when we have the initial time for this sub-segment

  last_sample=ii

! loop over the timesteps required in the resampled dataset  
  do i_time=1,nt_sub
  
    time=tstart+(i_time-1)*dt_sub
    
! loop over the input dataset and find the timesteps which bracket the required time
    do ii=last_sample,nt_pad-1
    
      t1=t(ii)
      t2=t(ii+1)
      
      if ( (time.GE.t1).AND.(time.LT.t2) ) then
! the required time is bracketed by timesteps ii and ii+1
        GOTO 7010
      end if
    
    end do
    
    if (time.EQ.t2) GOTO 7010    ! special case for the last sample
    
    write(*,*)'ERROR: Unable to find the timesteps bracketing time:',time
    STOP 1
      
7010 CONTINUE   ! Jump here when we have bracketed the required time

    value_1=ft_filtered(ii)
    value_2=ft_filtered(ii+1)
	   
    interpolated_value=value_1+( (value_2-value_1)/(t2-t1) )*(time-t1)
    ft_sub(i_time)=interpolated_value
	    
    last_sample=ii    
            
  end do ! next i_time


!4 WRITE SUB-SAMPLE DATA TO .WAV FILE

n_samples=nt_sub

real_max=-1e30
real_min=1e30
real_avg=0.0

do sample=1,n_samples
  
  real_max=max(ft_sub(sample),real_max)
  real_min=min(ft_sub(sample),real_min)
  real_avg=real_avg+ft_sub(sample)

end do

real_avg=real_avg/real(n_samples)

real_range=(real_max-real_min)

write(*,*)
write(*,*)'maximum waveform value=',real_max
write(*,*)'minimum waveform value=',real_min
write(*,*)'waveform value range  =',real_range
write(*,*)'average waveform value=',real_avg

write(*,*)
write(*,*)'subtract average value and convert to integer*2'
write(*,*)

! subtract the average value from the function

ft_sub(:)=ft_sub(:)-real_avg

! convert the dataset to integer*2

ALLOCATE( int_data(1:n_samples) )

do sample=1,n_samples
 
  int_data(sample)=NINT(ft_sub(sample)*65534.0/real_range) 
  
end do

! set the wav format information

 ChunkID='RIFF'
 Format='WAVE'
 Subchunk1ID='fmt '
 Subchunk1Size=16
 AudioFormat=1
 NumChannels=1
 SampleRate=44100
 ByteRate=88200
 BlockAlign=2
 BitsPerSample=16
 
 Subchunk2ID='data'
 Subchunk2Size=n_samples*2       ! two bytes per sample (integer*2 data)
 
 ChunkSize=Subchunk1Size+Subchunk2Size+20

! write wav format file

write(*,*)'Enter the name for the .wav file (without .wav extension)'
read(*,'(A)')filename

write(record_user_inputs_unit,'(A)')trim(filename)
write(post_process_info_unit,*)'.wav filename (without .wav extension):',trim(ipfilename)

opfilename=trim(filename)//'.wav'

open(unit=local_file_unit,file=opfilename,form='unformatted',access='stream')

write(local_file_unit)ChunkID
write(local_file_unit)ChunkSize
write(local_file_unit)Format
write(local_file_unit)Subchunk1ID
write(local_file_unit)Subchunk1Size
write(local_file_unit)AudioFormat
write(local_file_unit)NumChannels
write(local_file_unit)SampleRate
write(local_file_unit)ByteRate
write(local_file_unit)BlockAlign
write(local_file_unit)BitsPerSample
write(local_file_unit)Subchunk2ID
write(local_file_unit)Subchunk2Size

do i=1,n_samples
   write(local_file_unit) int_data(i)
enddo

close(unit=local_file_unit)

!11. FINISH OFF

DEALLOCATE( ft )
DEALLOCATE( t )
DEALLOCATE( frequency )
DEALLOCATE( ft_filtered )
DEALLOCATE( ft_impulse )
DEALLOCATE( t_sub )
DEALLOCATE( ft_sub )

RETURN

9000 write(*,*)'Error reading the input datafile'
     STOP 1

END SUBROUTINE convert_data_to_sound
