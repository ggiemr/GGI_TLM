SUBROUTINE post_process_time_domain_data

! process time domain conducted emissions data 
USE constants
USE file_information

IMPLICIT NONE

integer :: nt,nt_pad
character(LEN=2) :: ch2
character        :: ch
character*256    :: ipfilename
character*256    :: filename
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

real*8 :: scale
real*8 :: Power_t_in

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
real*8    :: tstart,time
integer   :: ntstart,istart
real*8    :: window_spacing

real*8  :: t1,t2
integer :: last_sample,i_time
real*8  :: interpolated_value,value_1,value_2

! dataset padding stuff
integer   :: pad_factor,pad_factor2,zero_pad_factor
character :: pad_type

! time domain window stuff
character :: window_fn
real*8    :: wi

! frequency domain output with specified bandwidth stuff
real*8  :: fmin,fmax,fstep,f
integer :: nf
real*8  :: BW
character :: detector_type
logical :: BW_warning
logical :: freq_range_warning
logical :: centre_freq_found
character*256    :: opfilename
real*8  :: detector_half_width
integer :: hw_samples,ifsample
real*8  :: fs,df

integer :: n_detectors,detector

real*8  :: f_detector,V_detector,P_detector
real*8  :: Pmultiplier
integer :: if1,if2,last_if1

real*8,allocatable     :: frequency_sub(:)
real*8,allocatable     :: t_sub(:)
complex*16,allocatable :: ft_sub(:,:)
complex*16,allocatable :: ft_sub_temp(:)
complex*16,allocatable :: ft_sub_tavg(:)
complex*16,allocatable :: ft_sub_favg(:)

complex*16 :: Vdc

integer i,ii,isegment
  
logical		:: file_exists

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

write(*,*)'Enter the column number for voltage data'
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

!2. SCALE THE DATA (E.G. TO MICRO VOLTS) IF REQUIRED

write(*,*)
write(*,*)'Enter the scaling factor for the data (eg 1e6 to scale from V to uV)'
read(*,*)scale
write(record_user_inputs_unit,*)scale,'  # scale factor for data (eg 1E6 for V to micro V'
write(*,*)'Scale factor is:',scale

ft(:)=ft(:)*scale
Pmultiplier=1.0/(scale*scale)

! Calculate the sum of the squares of the samples i.e. a measure of the power in the signal
Power_t_in=0d0
do i=1,nt
  Power_t_in=Power_t_in+ft(i)**2
end do

Power_t_in=power_t_in*dt*Pmultiplier

write(*,*)'Time domain input power:',Power_t_in

!3. APPLY LOW PASS AND HIGH PASS FILTERS TO INPUT TIME DOMAIN DATA AS REQUIRED

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

write(*,*)
write(*,*)'Input function:'
CALL Power_calc_time(ft_filtered,nt_pad,dt,Pmultiplier)

! Fourier Transform the input funtion
  
CALL FFT(ft_impulse,nt_pad)
CALL FFT(ft_filtered,nt_pad)

CALL Power_calc_FFT(ft_filtered,nt_pad,dt,Pmultiplier)
  
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
 
write(*,*)
write(*,*)'Filtered function:'
CALL Power_calc_FFT(ft_filtered,nt_pad,dt,Pmultiplier)

CALL INVERSE_FFT(ft_impulse,nt_pad)
CALL INVERSE_FFT(ft_filtered,nt_pad)
  
if (write_temp_files) then
  filename='filter_impulse_response.dat'
  CALL write_FFT_data(t,ft_impulse,nt_pad,filename,local_file_unit)
  
  filename='ft_out.dat'
  CALL write_FFT_data(t,ft_filtered,nt_pad,filename,local_file_unit)
end if
  
CALL Power_calc_time(ft_filtered,nt_pad,dt,Pmultiplier)

!4 EXTRACT SUB-SAMPLES FROM THE INPUT DATASET

write(*,*)'Enter the number of sub-data segments to extract from the original dataset'
read(*,*)n_sub_segments
write(record_user_inputs_unit,*)n_sub_segments,'  # Number of sub-data segments of data to process'
write(*,*)'Number of sub-segments=',n_sub_segments

if (n_sub_segments.LT.1) then
  write(*,*)'ERROR: The number of sub-segments should be at least 1'
  STOP 1
end if

write(*,*)'Total number of input samples/ number of sub segments=',nt/n_sub_segments
write(*,*)'Period of input signal/ number of sub segments=',(t(nt)-t(1))/n_sub_segments

write(*,*)'Enter the period of each sub-segment'
read(*,*)t_sub_period
write(record_user_inputs_unit,*)t_sub_period,'  # period of each sub-segment'
write(*,*)'period of sub-segment=',t_sub_period

write(*,*)'Enter the number of samples in each sub-segment (this should be a power of 2)'
read(*,*)nt_sub
write(record_user_inputs_unit,*)nt_sub,'  # number of samples in each sub-segment (this should be a power of 2)'
write(*,*)'number of samples in each sub-segment=',nt_sub

! work out the power of two greater than or equal to nt_sub
nt_sub_pad2=2
do while (nt_sub_pad2.LT.nt_sub)
  nt_sub_pad2=nt_sub_pad2*2
end do
if (nt_sub_pad2.NE.nt_sub) then
  write(*,*)'WARNING: The number of sub-samples is not a power of 2. The next power of 2 is:',nt_sub_pad2
end if
nt_sub2=nt_sub
nt_sub=nt_sub_pad2   ! make a power of 2

write(*,*)'Enter the time period over which the samples are to be distributed or enter 0 for continuous subsets'
write(*,*)'tmax=',tmax
read(*,*)sub_sample_period
write(record_user_inputs_unit,*)sub_sample_period,  &
 '  # time period over which the samples are to be distributed or enter 0 for continuous subsets'
write(*,*)'sub sample period=',sub_sample_period

if (sub_sample_period.GT.tmax) then
  write(*,*)'ERROR: sub sample period should be less than tmax'
  STOP 1
end if

write(*,*)'Enter the padding factor (set to 1 for no padding, pad factor should be a power of two)'
read(*,*)pad_factor
write(record_user_inputs_unit,*)pad_factor,'  # pad factor (set to 1 for no padding, pad factor should be a power of two)'
write(*,*)'padding factor=',pad_factor

if (pad_factor.LT.1) then
  write(*,*)'ERROR: The padding factor should be at least 1'
  STOP 1
end if 
  
pad_factor2=1
do while (pad_factor2.LT.pad_factor)
  pad_factor2=pad_factor2*2
end do
if (pad_factor2.NE.pad_factor) then
  write(*,*)'ERROR: The pad_factor is not a power of 2. The next power of 2 is:',pad_factor2
  STOP 1
end if

write(*,*)'Enter z to pad with zeros or p to create a periodic continuation of the sub-signal'
read(*,'(A)'),pad_type
if ( (pad_type.NE.'z').AND.(pad_type.NE.'p') ) then
  write(*,*)'Pad type should be z or p'
  STOP 1
end if
write(record_user_inputs_unit,'(A,A)')pad_type,' # z to pad with zeros or p to create a periodic continuation of the sub-signal'
write(*,*)'pad_type=',pad_type

if (pad_type.EQ.'z')then
  zero_pad_factor=pad_factor
else
  zero_pad_factor=1
end if

write(*,*)'Enter the number of time domain sub-segments to plot'
read(*,*)n_sub_plot
write(record_user_inputs_unit,*)n_sub_plot,'  # number of sub-segments to plot (max 9)'
write(*,*)'Number of sub-segments to plot=',n_sub_plot

write(*,*)"Enter the type of time domain window function to apply; 'r'=rectangular, 'h'=Hann'"
read(*,'(A)')window_fn
write(record_user_inputs_unit,'(A,A)')window_fn,'  # time domain window function (r=rectangular, h=Hann)'
write(*,*)'window function=',window_fn

if ( (window_fn.NE.'r').AND.(window_fn.NE.'R').AND.(window_fn.NE.'h').AND.(window_fn.NE.'H') ) then
  write(*,*)'Unknown window function type:',window_fn
else if (window_fn.EQ.'r') then
  window_fn='R'
else if (window_fn.EQ.'h') then
  window_fn='H'
end if

nt_sub_pad=nt_sub*pad_factor

write(*,*)'Period of each sub-sample,    t_sub_period=',t_sub_period
write(*,*)'Number of time samples (increased to power of 2),    nt_sub=',nt_sub
write(*,*)'Number of time samples including zero padding, nt_sub_pad=',nt_sub_pad

if ( (n_sub_plot.LT.0).OR.(n_sub_plot.GT.9) ) then
  write(*,*)'Number of sub-segment datasets to plot should be in the range 0 to 9'
  STOP 1
end if

if ( n_sub_plot.GT.n_sub_segments ) then
  write(*,*)'Number of sub-segment datasets to plot exceeds the number of sub-segments'
  STOP 1
end if

if (sub_sample_period.NE.0d0) then
  window_spacing=sub_sample_period/dble(n_sub_segments)
else
  window_spacing=t_sub_period
end if

dt_sub=t_sub_period/dble(nt_sub)       ! timestep for the sub-sampled segments

write(*,*)'Timestep for sub-sampled segments =',dt_sub

! Checks for consistency

if (sub_sample_period.GT.nt*dt) then
  write(*,*)'ERROR: the sub-sampling goes beyond the time period of the input data'
  write(*,*)'Sub_sampling period=',sub_sample_period
  write(*,*)'Input data period  =',nt*dt
  STOP 1
end if

ALLOCATE( frequency_sub(1:nt_sub_pad) )
CALL FFT_set_frequencies(frequency_sub,nt_sub_pad,dt_sub)

! set the time axis data for sub-segments incuding zero padding
ALLOCATE( t_sub(1:nt_sub_pad) ) 
do i=1,nt_sub_pad
  t_sub(i)=dble(i-1)*dt_sub
end do

ALLOCATE( ft_sub(1:nt_sub_pad,1:n_sub_segments) )
ft_sub(1:nt_sub_pad,1:n_sub_segments)=(0d0,0d0)

write(*,*)
write(*,*)'extract the sub-segments from the original filtered dataset'

do isegment=1,n_sub_segments
  
  tstart=(isegment-1)*window_spacing
  
  write(*,*)
  write(*,*)'Sub segment             : ',isegment
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
    ft_sub(i_time,isegment)=interpolated_value
	    
    last_sample=ii    
            
  end do ! next i_time

  if (pad_type.EQ.'p') then
! Add a periodic extension of the dataset if required

    write(*,*)'Add periodic extension to sub-sample',isegment
    write(*,*)'copy samples',1,' to',nt_sub,' to period(s)'
    do i=2,pad_factor
      write(*,*)'            ',1+(i-1)*nt_sub,' to',nt_sub+(i-1)*nt_sub
    end do

! loop over the timesteps of the initial signal
    do i_time=1,nt_sub
  
      do i=2,pad_factor

        ii=i_time+(i-1)*nt_sub
!        write(*,*)'copy from',i_time,' to ',ii,'ft=',ft_sub(i_time,isegment)
        ft_sub(ii,isegment)=ft_sub(i_time,isegment)

      end do  ! next padding period
    
    end do ! next time sample of initial signal
  
  end if  ! Add periodic extension

end do    ! next sub_segment

!4. APPLY A WINDOW FUNCTION
! Note: no action reqiuired for rectangular window
if (window_fn.EQ.'H') then

write(*,*)'Apply Hann window'

do isegment=1,n_sub_segments

  do i=1,nt_sub_pad  
  
    wi=0.5d0*( 1d0-cos(2d0*pi*(dble(i-1))/dble(nt_sub_pad -1)) )
    ft_sub(i,isegment)=ft_sub(i,isegment)*wi
    
  end do

end do
  
end if

! allocate further arrays for processing and averaging purposes

ALLOCATE( ft_sub_temp(1:nt_sub_pad) )
ALLOCATE( ft_sub_tavg(1:nt_sub_pad) )
ALLOCATE( ft_sub_favg(1:nt_sub_pad) )

ft_sub_temp(1:nt_sub_pad)=(0d0,0d0)
ft_sub_tavg(1:nt_sub_pad)=(0d0,0d0)
ft_sub_favg(1:nt_sub_pad)=(0d0,0d0)

!5. AVERAGE IN TIME

do isegment=1,n_sub_segments

  do i=1,nt_sub_pad
  
    ft_sub_tavg(i)=ft_sub_tavg(i)+ft_sub(i,isegment)
    
  end do

end do

ft_sub_tavg(:)=ft_sub_tavg(:)/dble(n_sub_segments)

!6. CALCULATE TIME DOMAIN POWER

write(*,*)
write(*,*)'Time domain sub-samples:'

do isegment=1,n_sub_segments

  ft_sub_temp(:)=0d0
  ft_sub_temp(1:nt_sub_pad)=ft_sub(1:nt_sub_pad,isegment)
  CALL Power_calc_time(ft_sub_temp,nt_sub_pad,dt_sub,Pmultiplier)

end do

write(*,*)
write(*,*)'Time domain averaged sub-samples:'
CALL Power_calc_time(ft_sub_tavg,nt_sub_pad,dt_sub,Pmultiplier)

!7. Remove d.c. if required

write(*,*)'Do you want to subtract the d.c. voltage (y/n)'
read(*,'(A)')ch
write(record_user_inputs_unit,'(A,A)')ch,'  # subtract d.c. (y or n)'
write(*,*)'Remove d.c.:',ch

if ((ch.eq.'y').OR.(ch.EQ.'Y')) then

  write(*,*)'Subtracting the d.c. voltage'

  do isegment=1,n_sub_segments

    Vdc=0d0
    do i=1,nt_sub_pad
      Vdc=Vdc+ft_sub(i,isegment)
    end do
    
    Vdc=Vdc/nt_sub_pad
    
    write(*,*)'Vdc=',Vdc
    
    do i=1,nt_sub_pad
      ft_sub(i,isegment)=ft_sub(i,isegment)-Vdc
    end do

! Check that the d.c. removal has worked OK
    Vdc=0d0
    do i=1,nt_sub_pad
      Vdc=Vdc+ft_sub(i,isegment)
    end do
    
    Vdc=Vdc/nt_sub_pad
    
    write(*,*)'Check, final Vdc=',Vdc

  end do
  
end if

if (write_temp_files) then
! Plot time domain datasets
  do isegment=1,n_sub_plot

    write(ch,'(I1)')isegment
    filename='sub_segment_time_'//ch//'.dat'
    write(*,*)'Writing time domain dataset:',trim(filename)
    ft_sub_temp(1:nt_sub_pad)=ft_sub(1:nt_sub_pad,isegment)
    CALL write_FFT_data(t_sub,ft_sub_temp,nt_sub_pad,filename,local_file_unit)
  
  end do

!  filename='sub_segment_time_average.dat'
!  CALL write_FFT_data(t_sub,ft_sub_tavg,nt_sub_pad,filename,local_file_unit)
end if

!8. FFT

write(*,*)

! FFT each of the sub-segment data sets
do isegment=1,n_sub_segments

  ft_sub_temp(:)=(0d0,0d0)
  ft_sub_temp(1:nt_sub_pad)=ft_sub(1:nt_sub_pad,isegment)

  write(*,*)'FFT sub-segment',isegment,' number of samples=',nt_sub_pad

  CALL FFT(ft_sub_temp,nt_sub_pad)

  ft_sub(1:nt_sub_pad,isegment)=ft_sub_temp(1:nt_sub_pad)
  
end do

! FFT the time domain average dataset
CALL FFT(ft_sub_tavg,nt_sub_pad)

!9. AVERAGE THE MAGNITUDES IN FREQUENCY DOMAIN

ft_sub_favg(:)=(0d0,0d0)
do isegment=1,n_sub_segments

  do i=1,nt_sub_pad
  
    ft_sub_favg(i)=ft_sub_favg(i)+abs(ft_sub(i,isegment))
    
  end do

end do
ft_sub_favg(:)=ft_sub_favg(:)/dble(n_sub_segments)

if (write_temp_files) then

! Plot frequency domain datasets
  do i=1,n_sub_plot

    write(ch,'(I1)')i
    filename='sub_segment_freq_'//ch//'.dat'
    ft_sub_temp(:)=(0d0,0d0)
    ft_sub_temp(1:nt_sub_pad)=ft_sub(1:nt_sub_pad,i)
    write(*,*)'Writing time domain dataset:',trim(filename)
    CALL write_single_side_FFT_data(frequency_sub,ft_sub_temp,nt_sub_pad,zero_pad_factor,filename,local_file_unit)
  
  end do
end if

!if (write_temp_files) then
!  filename='sub_segment_time_average_freq.dat'
!  CALL write_single_side_FFT_data(frequency_sub,ft_sub_tavg,nt_sub_pad,zero_pad_factor,filename,local_file_unit)
!end if

filename='sub_segment_freq_average.dat'

write(*,*)'Enter the filename for the raw FFT output data'
read(*,'(A)')opfilename
write(record_user_inputs_unit,'(A)')trim(opfilename)
CALL write_single_side_FFT_data(frequency_sub,ft_sub_favg,nt_sub_pad,zero_pad_factor,opfilename,local_file_unit)

!10. CALCULATE FREQUENCY DOMAIN POWER 

write(*,*)
write(*,*)'Frequency domain power calculations:'

write(*,*)
write(*,*)'Frequency domain averaged magnitude over sub-samples:'
CALL Power_calc_FFT(ft_sub_favg,nt_sub_pad,dt_sub,Pmultiplier)

write(*,*)
write(*,*)'Time domain average over sub-samples:'
CALL Power_calc_FFT(ft_sub_tavg,nt_sub_pad,dt_sub,Pmultiplier)

!11. Output frequency domain data over a specified frequency ranges with given bandwidth detectors

write(*,*)
write(*,*)'Enter the number of frequency domain bands with specific detector bandwidths to output'
read(*,*)n_detectors
write(record_user_inputs_unit,*)n_detectors,' # number of frequency domain bands with specific detector bandwidths to output'
write(*,*)'Number of detectors:',n_detectors

write(*,*)'Enter the output filename'
read(*,'(A)')opfilename
write(record_user_inputs_unit,'(A)')trim(opfilename)

open(unit=local_file_unit,file=trim(opfilename))

do detector=1,n_detectors

  write(*,*)'Bandwidth and detector specification number ',detector
  write(*,*)
  write(*,*)'Enter the minimum frequency to output (Hz)'
  read(*,*)fmin
  write(record_user_inputs_unit,*)fmin,' # minimum frequency for output'
  write(*,*)'minimum frequency for output=',fmin

  write(*,*)'Enter the maximum frequency to output (Hz)'
  read(*,*)fmax
  write(record_user_inputs_unit,*)fmax,' # maximum frequency for output'
  write(*,*)'maximum frequency for output=',fmax

  write(*,*)'Enter the number of frequencies to output '
  read(*,*)nf
  write(record_user_inputs_unit,*)nf,' # number of frequencies for output'
  write(*,*)'number of frequencies for output=',nf

  write(*,*)"Enter the detector filter frequency domain function ('gaussian' or 'rectangular')"
  read(*,'(A)')detector_type
  write(record_user_inputs_unit,'(A,A)')detector_type,' # detector frequency domain function (gaussian or rectangular)'
  write(*,*)'detector_type=',detector_type

  if ( (detector_type.NE.'r').AND.(detector_type.NE.'R').AND.(detector_type.NE.'g').AND.(detector_type.NE.'G') ) then
    write(*,*)'Unknown detector type:',detector_type
  else if (detector_type.EQ.'r') then
    detector_type='R'
  else if (detector_type.EQ.'g') then
    detector_type='G'
  end if

  write(*,*)'Enter the detector bandwidth (Hz)'
  read(*,*)BW
  write(record_user_inputs_unit,*)BW,' # detector bandwidth'
  write(*,*)'detector bandwidth',BW

  fstep=(fmax-fmin)/dble(nf-1)

  freq_range_warning=.FALSE.
  BW_warning=.FALSE.

  df=1d0/(nt_sub_pad*dt_sub)

  write(*,*)'Frequency step check',df,frequency_sub(2)-frequency_sub(1)

  if (BW.LT.df) then
    BW_warning=.TRUE.
  end if

  fstep=(fmax-fmin)/dble(nf-1)

  if (detector_type.EQ.'R') then
    detector_half_width=BW/2d0
  else if (detector_type.EQ.'G') then
    detector_half_width=BW*2d0
  end if

  write(*,*)'Detector bandwidth =',BW
  write(*,*)'Detector half width=',detector_half_width
  write(*,*)'FFT frequency step =',df

  hw_samples=NINT(detector_half_width/df)+1    ! ensure that this is at least 1

  write(*,*)'Number of samples in detector half-width=',hw_samples

  write(*,*)' Check on the evaluation of the detector filter function in the frequency domain'

  do i=-hw_samples,hw_samples
    df=dble(i)/(nt_sub_pad*dt_sub)
    CALL evaluate_filter_function(detector_type,BW,df,f_detector)         
    write(*,*)i,df,BW,f_detector
  end do

  write(*,*)

  last_if1=1
  do i=1,nf
  
    f=fmin+dble(i-1)*fstep      ! centre frequency for the detector
  
! find the frequency samples bracketing the centre frequency 
    centre_freq_found=.FALSE.
    do if1=last_if1,nt_sub_pad/2-1
  
      if2=if1+1
    
      if ( (frequency_sub(if1).LE.f).AND.(frequency_sub(if2).GE.f) ) then
! if1 and if2 bracket the centre frequency
        centre_freq_found=.TRUE.
        last_if1=if1
        EXIT   ! exit from the loop
      end if
  
    end do
  
    P_detector=0d0   
  
    if (centre_freq_found) then

! loop over the samples within the bandwidth of the detector
      do ii=-hw_samples,hw_samples
      
        ifsample=if1+ii     ! sample number
            
        if ( (ifsample.GE.1).AND.(ifsample.LE.(nt_sub_pad/2)) ) then
      
          fs=frequency_sub(ifsample)
          df=abs(f-fs)
        
! evaluate the detector filter function at this frequency offset        
          CALL evaluate_filter_function(detector_type,BW,df,f_detector)         
        
! add the contribution of this frequency to the detector power. Note factor of 2 to take account of negative frequencies
! There is a question about what to do about zero padding here...

          V_detector=abs(f_detector*ft_sub_favg(ifsample)*zero_pad_factor/nt_sub_pad)
        
! Note in the power calculation, divide by pad_factor as the signal is padded out with zeros which should
! be corrected for in the calculation

          P_detector=P_detector+2d0*((V_detector)**2)/zero_pad_factor       
               
        else
          freq_range_warning=.TRUE.    ! warn that the filter goes out of the range of the frequency domain data
        end if
    
      end do
        
      if (P_detector.GT.0d0) then
        write(local_file_unit,'(4ES16.6)')f,P_detector,10d0*log10(P_detector/50d0)+30d0,10d0*log10(P_detector)   ! note conversion to dBm for 50 ohm load in col 3.
      else
        write(local_file_unit,'(4ES16.6)')f,P_detector,-200d0,-200d0                       ! zero power so write very small dBm value
      end if
    
    else
  
      freq_range_warning=.TRUE.
      write(*,*)'Centre frequeny out of range',f 
    
    end if

  end do  ! next frequency
  
  write(local_file_unit,*)
  write(local_file_unit,*)

end do  ! next detector

close(unit=local_file_unit)

if (BW_warning) then
  write(*,*)'***************************************************************************'
  write(*,*)'WARNING: the detector bandwidth is less than the frequency step in the data'
  write(*,*)'***************************************************************************'   
end if

if (freq_range_warning) then
  write(*,*)'***************************************************************************'
  write(*,*)'WARNING: the frequency sweep with this detector bandwidth goes out of range'
  write(*,*)'***************************************************************************'   
end if

!11. FINISH OFF

DEALLOCATE( ft )
DEALLOCATE( t )
DEALLOCATE( frequency )
DEALLOCATE( ft_filtered )
DEALLOCATE( ft_impulse )
DEALLOCATE( frequency_sub )
DEALLOCATE( t_sub )
DEALLOCATE( ft_sub )
DEALLOCATE( ft_sub_temp )
DEALLOCATE( ft_sub_tavg )
DEALLOCATE( ft_sub_favg )

RETURN

9000 write(*,*)'Error reading the input datafile'
     STOP 1

END SUBROUTINE post_process_time_domain_data
!
! ___________________________________________________________________________________
!
!
SUBROUTINE evaluate_filter_function(detector_type,BW,df,f_detector)     

IMPLICIT NONE

character :: detector_type
real*8    :: BW
real*8    :: df
real*8    :: f_detector

! local_variables

real*8 :: k

! START

if (detector_type.EQ.'R') then

  if (abs(df).LE.BW/2d0) then
    f_detector=1d0
  else
    f_detector=0d0
  end if
  RETURN

else if (detector_type.EQ.'G') then

  k=BW/(1.1774d0)
  f_detector=exp(-df*df/(k*k))
  RETURN

else

  write(*,*)'ERROR in evaluate_filter_function: unknown detector type:',detector_type
  STOP

end if


END SUBROUTINE evaluate_filter_function

!
! These subroutines come from the GGI_TLM project
!
! SUBROUTINE FFT_set frequencies
!
! NAME
!    FFTT_set frequencies
!
! DESCRIPTION
!     calculate the frequencies output by the FFT subroutine
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/11/2013 CJS
!
!
SUBROUTINE FFT_set_frequencies(frequency,n_frequencies,dt)

IMPLICIT NONE
  integer	:: n_frequencies
  real*8        :: frequency(1:n_frequencies)
  real*8        :: dt

! local variables

  real*8	:: fmin,fmax,fstep,fnyquist
  integer	:: frequency_loop
  
! START

  fmin=0d0
  fstep=1d0/(n_frequencies*dt)
  fmax=fstep*(n_frequencies-1)
  fnyquist=fstep*(n_frequencies/2d0)
  
  write(*,*)'Fmin    =',fmin,' Hz'
  write(*,*)'Fmax    =',fmax,' Hz'
  write(*,*)'FNyquist=',fnyquist,' Hz'
  write(*,*)'n_frequencies=',n_frequencies
  
  do frequency_loop=1,n_frequencies
      
    frequency(frequency_loop)=fmin+(frequency_loop-1)*fstep
    if (frequency(frequency_loop).GT.fnyquist+fstep/2) then
      frequency(frequency_loop)=frequency(frequency_loop)-(2d0*fnyquist)
    end if
  
  end do
   
  RETURN
  
END SUBROUTINE FFT_set_frequencies
!
! SUBROUTINE write_FFT_time_data
!
! NAME
!    write_FFT_time_data
!
! DESCRIPTION
!     write time domain data from a complex dataset from the FFT process
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
SUBROUTINE write_FFT_time_data(xdata,ydata,n,name,op_unit)


IMPLICIT NONE

  integer	:: n
  real*8        :: xdata(n)
  complex*16	:: ydata(n)
  character*256 :: name
  integer :: op_unit

! local variables

  integer 	:: i

! START

  open(unit=op_unit,file=trim(name))
  
! look at the first and last xdata samples to decide if we have frequency or time domain data

  if(xdata(n).GT.xdata(1)) then

! we have time domain data and should plot in order
    do i=1,n
      write(op_unit,'(4ES16.6)')xdata(i),real(ydata(i)),imag(ydata(i)),abs(ydata(i))
    end do
  
  else
! we have frequency domain data 
    write(*,*)'ERROR in write_FFT_time_data. It looks like it has been called with frequency domain data'
    STOP 1
  
  end if
  
  close(unit=op_unit)
  
  RETURN
 
END SUBROUTINE write_FFT_time_data
!
! SUBROUTINE write_FFT_data
!
! NAME
!    write_FFT_data
!
! DESCRIPTION
!     write a complex dataset from the FFT process. This could be time or frequency domain data.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
SUBROUTINE write_FFT_data(xdata,ydata,n,name,op_unit)


IMPLICIT NONE

  integer	:: n
  real*8        :: xdata(n)
  complex*16	:: ydata(n)
  character*256 :: name
  integer :: op_unit

! local variables

  integer 	:: i,n2

! START

  open(unit=op_unit,file=trim(name))
  
! look at the first and last xdata samples to decide if we have frequency or time domain data

  if(xdata(n).GT.xdata(1)) then

! we have time domain data and should plot in order
    write(*,*)'CALLED write_FFT_data: Writing time domain data'
    do i=1,n
      write(op_unit,'(4ES16.6)')xdata(i),real(ydata(i)),imag(ydata(i)),abs(ydata(i))
    end do
  
  else
! we have frequency domain data so reorder the data when writing and scale data by 1/n
    write(*,*)'CALLED write_FFT_data: Writing frequency domain data'
    n2=n/2
    do i=n2+2,n
      write(op_unit,'(4ES16.6)')xdata(i),real(ydata(i)/n),imag(ydata(i)/n),abs(ydata(i)/n)
    end do
    do i=1,n2+1
      write(op_unit,'(4ES16.6)')xdata(i),real(ydata(i)/n),imag(ydata(i)/n),abs(ydata(i)/n)
    end do
  
  end if
  
  close(unit=op_unit)
  
  RETURN
 
END SUBROUTINE write_FFT_data
!
! SUBROUTINE write_single_side_FFT_data
!
! NAME
!    write_single_side_FFT_data
!
! DESCRIPTION
!     write a complex dataset from the FFT process but only for positive frequencies
!     The output includes power into a 50 ohm system expressed in dBm for comparison with Spectrum Analyser output
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
SUBROUTINE write_single_side_FFT_data(xdata,ydata,n,pad_factor,name,op_unit)


IMPLICIT NONE

  integer	:: n
  real*8        :: xdata(n)
  complex*16	:: ydata(n)
  integer	:: pad_factor
  character*256 :: name
  integer :: op_unit

! local variables

  integer 	:: i,n2
  real*8        :: Vdbm
  real*8        :: Vdb

! START

  open(unit=op_unit,file=trim(name))
  
! look at the first and last xdata samples to decide if we have frequency or time domain data

  if(xdata(n).GT.xdata(1)) then

! we have time domain data 
    write(*,*)'CALLED write_single_side_FFT_data: Writing time domain data'
    write(*,*)'ERROR in write_single_side_FFT_data. It looks like it has been called with time domain data'
    STOP 1
  
  else
! we have frequency domain data so reorder the data when writing and scale data by 1/n
    write(*,*)'CALLED write_single_side_FFT_data: Writing frequency domain data'
    n2=n/2
    do i=1,n2+1
! power into 50ohm load in mW; note add 3dB so we can compare with single sided data
! Note also the use of the zero pad factor to scale 

      if (ydata(i).NE.0d0) then

        Vdbm=20d0*log10(abs(ydata(i)*pad_factor/n)/sqrt(50d0))+30d0+3.01d0      
        Vdb=20d0*log10(abs(ydata(i)*pad_factor/n))+3.01d0  
      
      else
      
        Vdbm=-200d0
        Vdb=-200d0
      
      end if
      
      write(op_unit,'(6ES16.6)')xdata(i),real(ydata(i)*pad_factor/n),imag(ydata(i)*pad_factor/n), &
                                     abs(ydata(i)*pad_factor/n),Vdbm,Vdb
    end do

  end if
  
  close(unit=op_unit)
  
  RETURN
 
END SUBROUTINE write_single_side_FFT_data
!
! SUBROUTINE butterworth_LPF
!
! NAME
!    
!
! DESCRIPTION
!     apply Butterworth low pass filter to frequency domain data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
SUBROUTINE butterworth_LPF(freq,v,n,order,fc,dt)


IMPLICIT NONE

  integer	:: n
  real*8        :: freq(n)
  complex*16	:: v(n)
  integer       :: order
  real*8        :: fc
  real*8        :: dt

! START

! local variables

  integer 	:: i,n2
  real*8        :: w,wc
  real*8        :: wt,wct
  complex*16    :: H,s
  
  complex*16,parameter :: j=(0d0,1d0)
  real*8,parameter :: pi=3.1415926535d0

! START

  wc=fc*2d0*pi
! apply the frequency warping function to wc
  wct=(2d0/dt)*tan(wc*dt/2d0)
    
  n2=n/2

  do i=1,n
    
    if ( (i.NE.1).AND.(i.NE.n2+1) )  then   ! these are not zero or maximum frequency terms
      
      w=2d0*pi*freq(i)
! apply the frequency warping function to w
      wt=(2d0/dt)*tan(w*dt/2d0)
      
      s=j*wt/wct         ! normalised Laplace variable
    
      if(order.EQ.1) then
    
        H=1d0/(s+1d0)
      
      else if(order.EQ.2) then
    
        H=1d0/(s*s+1.414213562d0*s+1)
      
      else if(order.EQ.3) then
    
        H=1d0/((s+1)*(s*s+s+1))
      
      else if(order.EQ.4) then
    
        H=1d0/((s*s+0.7654d0*s+1)*(s*s+1.8478d0*s+1))
      
      else
        write(*,*)'No implementation for Butterworth low pass filter of order:',order
        STOP 1
      end if
      
    else if (i.EQ.n2+1) then    ! set H as the value as s-> infinity
    
      H=(0d0,0d0)
    
    else      ! set H as the value as s-> zero
    
      H=(1d0,0d0)
    
    end if
    
    v(i)=v(i)*H
    
  end do
  
  
  RETURN
 
END SUBROUTINE butterworth_LPF
!
! SUBROUTINE butterworth_HPF
!
! NAME
!    
!
! DESCRIPTION
!     apply Butterworth high pass filter to frequency domain data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
SUBROUTINE butterworth_HPF(freq,v,n,order,fc,dt)


IMPLICIT NONE

  integer	:: n
  real*8        :: freq(n)
  complex*16	:: v(n)
  integer       :: order
  real*8        :: fc
  real*8        :: dt

! local variables

  integer 	:: i,n2
  real*8        :: w,wc
  real*8        :: wt,wct
  complex*16    :: H,s
  
  complex*16,parameter :: j=(0d0,1d0)
  real*8,parameter :: pi=3.1415926535d0

! START

  wc=fc*2d0*pi
! apply the frequency warping function to wc
  wct=(2d0/dt)*tan(wc*dt/2d0)
  
  n2=n/2

  do i=1,n
    
    if ( (i.NE.1).AND.(i.NE.n2+1) )  then   ! these are not zero or maximum frequency terms
      
      w=2d0*pi*freq(i)
! apply the frequency warping function to w
      wt=(2d0/dt)*tan(w*dt/2d0)
    
      s=j*wt/wct         ! normalised Laplace variable
      s=1d0/s         ! transformation which turns LPF filter functions into HPF filter functions
    
      if(order.EQ.1) then
    
        H=1d0/(s+1d0)
      
      else if(order.EQ.2) then
    
        H=1d0/(s*s+1.414213562d0*s+1)
      
      else if(order.EQ.3) then
    
        H=1d0/((s+1)*(s*s+s+1))
      
      else if(order.EQ.4) then
    
        H=1d0/((s*s+0.7654d0*s+1)*(s*s+1.8478d0*s+1))
      
      else
        write(*,*)'No implementation for Butterworth low pass filter of order:',order
        STOP 1
      end if
      
    
    else if (i.EQ.n2+1) then    ! set H as the value as s-> infinity
    
      H=(1d0,0d0)
    
    else      ! set H as the value as s-> zero
    
      H=(0d0,0d0)
    
    end if
    
    v(i)=v(i)*H
    
  end do
  
  
  RETURN
 
END SUBROUTINE butterworth_HPF
!
! SUBROUTINE Power_calc_FFT
!
! NAME
!    
!
! DESCRIPTION
!     calculate the average power from frequency domain data
!     this is calculated as the integral of v^2 over frequency
!
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
SUBROUTINE Power_calc_FFT(v,n,dt,Pmultiplier)


IMPLICIT NONE

  integer	:: n
  complex*16	:: v(n)
  real*8        :: dt,Pmultiplier

! local variables

  integer 	:: i
  real*8        :: Power
  
  complex*16,parameter :: j=(0d0,1d0)
  real*8,parameter :: pi=3.1415926535d0

! START

  Power=0d0
    
  do i=1,n
    
    Power=Power+v(i)*conjg(v(i))
    
  end do

  Power=Power*Pmultiplier/(dble(n)*dble(n))

  write(*,*)'Power calculated from Frequency domain data:'
  write(*,*)Power,' ',10d0*log10(Power),'dB'
  
  RETURN
 
END SUBROUTINE Power_calc_FFT
!
! SUBROUTINE Power_calc_time
!
! NAME
!    
!
! DESCRIPTION
!     calculate the average power from time domain data
!     this is calculated as the average of v^2 in time
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
SUBROUTINE Power_calc_time(v,n,dt,Pmultiplier)


IMPLICIT NONE

  integer	:: n
  complex*16	:: v(n)
  real*8        :: dt,Pmultiplier

! local variables

  integer 	:: i
  real*8        :: Power
  
  complex*16,parameter :: j=(0d0,1d0)
  real*8,parameter :: pi=3.1415926535d0

! START

  Power=0d0
  
  do i=1,n
    
    Power=Power+dble(V(i))*conjg(v(i))*dt*Pmultiplier
    
  end do
  
  Power=Power/(n*dt)

  write(*,*)'Power calculated from Time domain data:'
  write(*,*)Power,' ',10d0*log10(Power),'dB'

    
  RETURN
 
END SUBROUTINE Power_calc_time
