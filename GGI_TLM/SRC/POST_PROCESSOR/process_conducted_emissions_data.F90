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
integer :: window_spacing
integer :: n_sub_plot
integer :: zero_pad_factor
integer :: nt_sub_pad,nt_sub_pad2
real*8    :: sub_sample_period
real*8    :: dt_sub

! time domain window stuff
character :: window_fn
real*8    :: wi

! frequency domain output with specified bandwidth stuff
real*8 :: fmin,fmax,fstep,f
integer :: nf
real*8 :: BW
character :: detector_type
logical :: BW_warning
logical :: freq_range_warning
logical :: centre_freq_found
character*256    :: opfilename
real*8  :: detector_half_width
integer :: hw_samples,ifsample
real*8  :: fs,df
real*8  :: f_detector,V_detector,P_detector
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

logical,parameter :: write_temp_files=.FALSE.

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

write(*,*)'Number of samples read, nt=',nt
write(*,*)'Number of samples with zero padding, nt_pad=',nt_pad

CLOSE(unit=local_file_unit)
DEALLOCATE( data_line )

dt=t(2)-t(1)

! fill the time axis data including the padded values
do i=1,nt_pad
  t(i)=dble(i-1)*dt
end do

!2. SCALE THE DATA (E.G. TO MICRO VOLTS) IF REQUIRED

write(*,*)'Enter the scaling factor for the data (eg 1e6 to scale from V to uV)'
read(*,*)scale
write(record_user_inputs_unit,*)scale,'  # scale factor for data (eg 1E6 for V to micro V'

ft(:)=ft(:)*scale

! Calculate the sum of the squares of the samples i.e. a measure of the power in the signal
Power_t_in=0d0
do i=1,nt
  Power_t_in=Power_t_in+ft(i)**2
end do

Power_t_in=power_t_in*dt

write(*,*)'Time domain input power:',Power_t_in,' W'

!3. APPLY LOW PASS AND HIGH PASS FILTERS TO INPUT TIME DOMAIN DATA AS REQUIRED

! Work out the frequency for each sample
ALLOCATE( frequency(1:nt_pad) )
CALL FFT_set_frequencies(frequency,nt_pad,dt)

ALLOCATE( ft_filtered(1:nt_pad) )

write(*,*)'Enter the number of filters to apply to the input time domain data'
read(*,*)n_filter
write(record_user_inputs_unit,*)n_filter,'  # number of filters to apply'

ft_filtered(1:nt_pad)=dcmplx(ft(1:nt_pad))

! allocate data for the impulse response
ALLOCATE( ft_impulse(1:nt_pad) )
ft_impulse(1:nt_pad)=(0d0,0d0)
ft_impulse(1)=(1d0,0d0)

write(*,*)
write(*,*)'Input function:'
CALL Power_calc_time(ft_filtered,nt_pad,dt)

! Fourier Transform the input funtion
  
CALL FFT(ft_impulse,nt_pad)
CALL FFT(ft_filtered,nt_pad)

CALL Power_calc_FFT(ft_filtered,nt_pad,dt)
  
if (write_temp_files) then
  filename='ft_freq.dat'
  CALL write_FFT_data(frequency,ft_filtered,nt_pad,filename,local_file_unit)
end if

do filter=1,n_filter

  write(*,*)'Enter the filter type (LPF or HPF)'
  read(*,*)ch
  write(record_user_inputs_unit,'(A,A)')ch,'  # filter type (LPF or HPF)'

  write(*,*)'Enter the filter order (integer)'
  read(*,*)order
  write(record_user_inputs_unit,*)order,'  # filter order'
  
  write(*,*)'Enter the filter 3dB frequency (Hz)'
  read(*,*)fc
  write(record_user_inputs_unit,*)fc,'  # filter cutoff frequency'
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
CALL Power_calc_FFT(ft_filtered,nt_pad,dt)

CALL INVERSE_FFT(ft_impulse,nt_pad)
CALL INVERSE_FFT(ft_filtered,nt_pad)
  
if (write_temp_files) then
  filename='filter_impulse_response.dat'
  CALL write_FFT_data(t,ft_impulse,nt_pad,filename,local_file_unit)
  
  filename='ft_out.dat'
  CALL write_FFT_data(t,ft_filtered,nt_pad,filename,local_file_unit)
end if
  
CALL Power_calc_time(ft_filtered,nt_pad,dt)

!4 EXTRACT SUB-SAMPLES FROM THE INPUT DATASET

write(*,*)'Enter the number of sub-data segments to extract from the original dataset'
read(*,*)n_sub_segments
write(record_user_inputs_unit,*)n_sub_segments,'  # Number of sub-data segments of data to process'

write(*,*)'Enter the number of samples in each sub-segment (this should be a power of 2)'
read(*,*)nt_sub
write(record_user_inputs_unit,*)nt_sub,'  # number of samples in each sub-segment (this should be a power of 2)'

! work out the power of two greater than or equal to nt
nt_sub_pad2=2
do while (nt_sub_pad2.LT.nt_sub)
  nt_sub_pad2=nt_sub_pad2*2
end do
if (nt_sub_pad2.NE.nt_sub) then
  write(*,*)'WARNING: The number of sub-samples is not a power of 2. The next power of 2 is:',nt_sub_pad2
  write(*,*)'Zero padding will be required for the FFT'
end if
nt_sub2=nt_sub
nt_sub=nt_sub_pad2   ! make a power of 2

write(*,*)'Enter the step size for sampling the original dataset (allows lower sample rate)'
read(*,*)step_sub
write(record_user_inputs_unit,*)step_sub,'  # step size for sub-segment sampling. This allows sampling at a lower sample rate'

write(*,*)'Enter the time period over which the samples are to be distributed or enter 0 for continuous subsets'
read(*,*)sub_sample_period
write(record_user_inputs_unit,*)sub_sample_period,  &
 '  # time period over which the samples are to be distributed or enter 0 for continuous subsets'

write(*,*)'Enter the zero padding factor (should be a power of two)'
read(*,*)zero_pad_factor
write(record_user_inputs_unit,*)zero_pad_factor,'  # zero pad factor (should be a power of two)'

write(*,*)'Enter the number of time domain sub-segments to plot'
read(*,*)n_sub_plot
write(record_user_inputs_unit,*)n_sub_plot,'  # number of sub-segments to plot (max 9)'

write(*,*)"Enter the type of time domain window function to apply; 'r'=rectangular, 'h'=Hann'"
read(*,'(A)')window_fn
write(record_user_inputs_unit,'(A,A)')window_fn,'  # time domain window function (r=rectangular, h=Hann)'

if ( (window_fn.NE.'r').AND.(window_fn.NE.'R').AND.(window_fn.NE.'h').AND.(window_fn.NE.'H') ) then
  write(*,*)'Unknown window function type:',window_fn
else if (window_fn.EQ.'r') then
  window_fn='R'
else if (window_fn.EQ.'h') then
  window_fn='H'
end if

nt_sub_pad=nt_sub*zero_pad_factor

write(*,*)'Number of time samples in each sub-sample,         n_sub2=',nt_sub2
write(*,*)'Number of time samples increased to power of 2,    nt_sub=',nt_sub
write(*,*)'Number of time samples including zero padding, nt_sub_pad=',nt_sub_pad

if ( (n_sub_plot.LT.0).OR.(n_sub_plot.GT.9) ) then
  write(*,*)'Number of sub-segment datasets to plot should be in the range 0 to 9'
  STOP 1
end if

if ( n_sub_plot.GT.n_sub_segments ) then
  write(*,*)'Number of sub-segment datasets to plot exceeds the number of sub-segments'
  STOP 1
end if

if (window_spacing.NE.0d0) then
  window_spacing=NINT(sub_sample_period/dble(n_sub_segments))
else
  window_spacing=nt_sub
end if

dt_sub=dt*step_sub       ! timestep for the sub-sampled segments

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

write(*,*)'extract the sub-segments from the original filtered dataset'

do isegment=1,n_sub_segments

  write(*,*)'First sample of sub_segment is:',(isegment-1)*window_spacing

  do i=1,nt_sub2   ! note use the number of samples requested initially, not the number of samples with zero padding
  
    ii=(isegment-1)*window_spacing+1+(i-1)*step_sub
    if (ii.GT.nt) then
      write(*,*)'ERROR: sub-sampling beyond the input data range'
      write(*,*)'Data sample:',ii,' Total number of samples=',nt
      STOP 1
    end if
    ft_sub(i,isegment)=ft_filtered(ii)
    
  end do

end do

!4. APPLY A WINDOW FUNCTION
! Note: no action reqiuired for rectangular window
if (window_fn.EQ.'H') then

do isegment=1,n_sub_segments

  do i=1,nt_sub2  ! note use the number of samples requested initially, not the number of samples with zero padding
  
    wi=0.5d0*( 1d0-cos(2d0*pi*(dble(i-1))/dble(nt_sub2-1)) )
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

if (write_temp_files) then
! Plot time domain datasets
  do isegment=1,n_sub_plot

    write(ch,'(I1)')isegment
    filename='sub_segment_time_'//ch//'.dat'
    ft_sub_temp(1:nt_sub)=ft_sub(1:nt_sub,isegment)
    CALL write_FFT_data(t_sub,ft_sub_temp,nt_sub_pad,filename,local_file_unit)
  
  end do

  filename='sub_segment_time_average.dat'
  CALL write_FFT_data(t_sub,ft_sub_tavg,nt_sub_pad,filename,local_file_unit)
end if

!6. CALCULATE TIME DOMAIN POWER

write(*,*)
write(*,*)'Time domain sub-samples:'

do isegment=1,n_sub_segments

  ft_sub_temp(:)=0d0
  ft_sub_temp(1:nt_sub_pad)=ft_sub(1:nt_sub_pad,isegment)
  CALL Power_calc_time(ft_sub_temp,nt_sub_pad,dt_sub)

end do

write(*,*)
write(*,*)'Time domain averaged sub-samples:'
CALL Power_calc_time(ft_sub_tavg,nt_sub_pad,dt_sub)

!7. Remove d.c. if required

write(*,*)'Do you want to subtract the d.c. voltage (y/n)'
read(*,'(A)')ch
write(record_user_inputs_unit,'(A,A)')ch,'  # subtract d.c. (y or n)'

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

!8. FFT

! FFT each of the sub-segment data sets
do isegment=1,n_sub_segments

  ft_sub_temp(:)=(0d0,0d0)
  ft_sub_temp(1:nt_sub_pad)=ft_sub(1:nt_sub_pad,isegment)

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
    CALL write_single_side_FFT_data(frequency_sub,ft_sub_temp,nt_sub_pad,zero_pad_factor,filename,local_file_unit)
  
  end do
end if

if (write_temp_files) then
  filename='sub_segment_time_average_freq.dat'
  CALL write_single_side_FFT_data(frequency_sub,ft_sub_tavg,nt_sub_pad,zero_pad_factor,filename,local_file_unit)
end if

filename='sub_segment_freq_average.dat'
CALL write_single_side_FFT_data(frequency_sub,ft_sub_favg,nt_sub_pad,zero_pad_factor,filename,local_file_unit)

!10. CALCULATE FREQUENCY DOMAIN POWER 

write(*,*)
write(*,*)'Frequency domain power calculations:'

write(*,*)
write(*,*)'Frequency domain averaged magnitude over sub-samples:'
CALL Power_calc_FFT(ft_sub_favg,nt_sub_pad,dt_sub)

write(*,*)
write(*,*)'Time domain average over sub-samples:'
CALL Power_calc_FFT(ft_sub_tavg,nt_sub_pad,dt_sub)

!11. Output frequency domain data over a specified frequency range with a given bandwidth detector
write(*,*)
write(*,*)'Enter the minimum frequency to output (Hz)'
read(*,*)fmin
write(record_user_inputs_unit,*)fmin,' # minimum frequency for output'

write(*,*)'Enter the maximum frequency to output (Hz)'
read(*,*)fmax
write(record_user_inputs_unit,*)fmax,' # maximum frequency for output'

write(*,*)'Enter the number of frequencies to output '
read(*,*)nf
write(record_user_inputs_unit,*)nf,' # number of frequencies for output'

write(*,*)"Enter the detector filter frequency domain function ('gaussian' or 'rectangular')"
read(*,'(A)')detector_type
write(record_user_inputs_unit,'(A,A)')detector_type,' # detector frequency domain function (gaussian or rectangular)'

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

write(*,*)'Enter the output filename'
read(*,'(A)')opfilename
write(record_user_inputs_unit,'(A)')trim(opfilename)

fstep=(fmax-fmin)/dble(nf-1)

open(unit=local_file_unit,file=trim(opfilename))

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
        
! Note in the power calculation, divide by zero_pad_factor as the signal is padded out with zeros which should
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

end do

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
    do i=1,n
      write(op_unit,'(4ES16.6)')xdata(i),real(ydata(i)),imag(ydata(i)),abs(ydata(i))
    end do
  
  else
! we have frequency domain data so reorder the data when writing and scale data by 1/n
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
SUBROUTINE write_single_side_FFT_data(xdata,ydata,n,zero_pad_factor,name,op_unit)


IMPLICIT NONE

  integer	:: n
  real*8        :: xdata(n)
  complex*16	:: ydata(n)
  integer	:: zero_pad_factor
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
    write(*,*)'ERROR in write_single_side_FFT_data. It looks like it has been called with time domain data'
    STOP 1
  
  else
! we have frequency domain data so reorder the data when writing and scale data by 1/n
    n2=n/2
    do i=1,n2+1
! power into 50ohm load in mW; note add 3dB so we can compare with single sided data
! Note also the use of the zero pad factor to scale 

      Vdbm=20d0*log10(abs(ydata(i)*zero_pad_factor/n)/sqrt(50d0))+30d0+3.01d0      
      Vdb=20d0*log10(abs(ydata(i)*zero_pad_factor/n))+3.01d0  
      
      write(op_unit,'(6ES16.6)')xdata(i),real(ydata(i)*zero_pad_factor/n),imag(ydata(i)*zero_pad_factor/n), &
                                     abs(ydata(i)*zero_pad_factor/n),Vdbm,Vdb
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
SUBROUTINE Power_calc_FFT(v,n,dt)


IMPLICIT NONE

  integer	:: n
  complex*16	:: v(n)
  real*8        :: dt

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

  Power=Power/(dble(n)*dble(n))

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
SUBROUTINE Power_calc_time(v,n,dt)


IMPLICIT NONE

  integer	:: n
  complex*16	:: v(n)
  real*8        :: dt

! local variables

  integer 	:: i
  real*8        :: Power
  
  complex*16,parameter :: j=(0d0,1d0)
  real*8,parameter :: pi=3.1415926535d0

! START

  Power=0d0
  
  do i=1,n
    
    Power=Power+dble(V(i))*conjg(v(i))*dt
    
  end do
  
  Power=Power/(n*dt)

  write(*,*)'Power calculated from Time domain data:'
  write(*,*)Power,' ',10d0*log10(Power),'dB'

    
  RETURN
 
END SUBROUTINE Power_calc_time
