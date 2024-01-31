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
! SUBROUTINE EMI_receiver_model
!
! NAME  EMI_receiver_model
!    
!
! DESCRIPTION
!      EMI_receiver_model based on 
! L. Yang, S. Wang, H. Zhao, Y. Zhi, "Prediction and Analysis of EMI Spectrum Based
!      on the Operating principle of EMC Spectrum Analyzers," IEEE trans Power Electronics,
!      Vol 35, No 1, January 2020.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 31/01/2024 CJS
!
!
SUBROUTINE EMI_receiver_model

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: time_function_number=1
  integer	:: freq_function_number=1
  integer	:: EMI_function_number=2
 
  real*8  :: FSWEEP_fmin,FSWEEP_fmax,FSWEEP_fstep
  integer :: n_FSWEEP

  real*8  :: RBW,FILTER_c
  real*8  :: FFILTER_fmin,FFILTER_fmax,FFILTER_fstep
  integer :: FFILTER_nf

  real*8,allocatable     :: filter_response(:)
  complex*16,allocatable :: filtered_signal(:)
  real*8,allocatable     :: Amag(:),Aphase(:),Aw(:)
  integer                :: nfsample

  integer :: n_envelope
  real*8  :: tmax_envelope,dt_envelope,t
  real*8,allocatable :: envelope(:)

  real*8 :: stemp
  real*8 :: Aw0,CT,ST

  integer :: nf0,nfm,nfp

  real*8 :: Vpeak,Vrms,Vavg,Vquasi_peak,Vout

  integer :: i,ii,ij,ik
  integer :: op_period

  real*8  :: Vqp,Vqpmin,Vqpmax,LHS,RHS
  
  integer :: floop,tloop
  
  integer :: detector_type

  real*8 :: TD,TC
  
  real*8 :: tmax

  real*8	:: fmin,fmax,fstep,f
  integer	:: n_frequencies
    
  integer	:: frequency_loop
  integer	:: timestep
  
  real*8	:: dt
  real*8	:: w
  complex*16	:: integral
  
  character*3	:: freq_range_type
  
  integer,parameter :: detector_type_peak=1
  integer,parameter :: detector_type_rms=2
  integer,parameter :: detector_type_quasi_peak=3
  integer,parameter :: detector_type_avg=4
  
! START

  write(*,*)
  write(*,*)'Apply EMI receiver model to time domain data'
  write(*,*)'NOTE: This process is experimental at the moment and requires proper testing and validation...'
  write(*,*)

  write(post_process_info_unit,*)
  write(post_process_info_unit,*)'Apply EMI receiver model to time domain data'
  write(post_process_info_unit,*)'NOTE: This process is experimental at the moment and requires proper testing and validation...'
  write(post_process_info_unit,*)

  n_functions_of_time=1
  n_functions_of_frequency=2
  
  CALL Allocate_post_data()
  
  CALL read_Time_Domain_Data(time_function_number) ! read first and only function of time
  
  dt=function_of_time(time_function_number)%time(2)-function_of_time(time_function_number)%time(1)
  write(*,*)'Timestep=',dt
  write(post_process_info_unit,*)'Timestep=',dt

! setup the frequency list

100 CONTINUE

  freq_range_type='lin'

! Specify the EMI receiver model frequency sweep

  write(*,*)'Enter minimum frequency for EMI receiver sweep, fmin'
  read(*,*)FSWEEP_fmin
  write(*,*)'Enter maximum frequency for EMI receiver sweep, fmax'
  read(*,*)FSWEEP_fmax
  write(*,*)'Enter the number of frequencies for EMI receiver sweep'
  read(*,*)n_FSWEEP
  
  write(record_user_inputs_unit,'(E16.7,A)')FSWEEP_fmin,' EMI frequency sweep fmin'
  write(record_user_inputs_unit,'(E16.7,A)')FSWEEP_fmax,' EMI frequency sweep fmax'
  write(record_user_inputs_unit,'(I16,A)')n_FSWEEP,' EMI frequency sweep n_frequencies'
  
  write(post_process_info_unit,*)'EMI frequency sweep Fmin=',FSWEEP_fmin,' Hz'
  write(post_process_info_unit,*)'EMI frequency sweep Fmax=',FSWEEP_fmax,' Hz'
  write(post_process_info_unit,*)'EMI frequency sweep nf=',n_FSWEEP
    
  if (n_frequencies.ne.1) then
    FSWEEP_fstep=(FSWEEP_fmax-FSWEEP_fmin)/dble(n_FSWEEP-1)
  else
    FSWEEP_fstep=0d0
  end if
  
! Specify the EMI receiver bandwidth

  write(*,*)'Enter EMI receiver resolution bandwidth, RBW'
  read(*,*)RBW
  
  write(record_user_inputs_unit,'(E16.7,A)')RBW,' EMI receiver resolution bandwidth, RBW'
  write(post_process_info_unit,*)'EMI receiver resolution bandwidth, RBW=',RBW,' Hz'
 
  write(*,*)'Enter EMI receiver detector type:'
  write(*,*)'1: Peak detector'
  write(*,*)'2: RMS detector'
  write(*,*)'3: Quasi-Peak detector'
  write(*,*)'4: Average detector'
  read(*,*)detector_type
  
  if ( (detector_type.LT.1).OR.(detector_type.GT.4) ) then
    write(*,*)'ERROR: the detector type should be in the range 1 to 4'
    STOP 1
  end if
  
  write(record_user_inputs_unit,'(I4,A)')detector_type,' EMI detector type'
  write(post_process_info_unit,*)'EMI detector_type=',detector_type
  
  if (detector_type.EQ.detector_type_quasi_peak) then
  
    write(*,*)'Enter quasi-peak discharge time constant, TD:'
    read(*,*)TD
    write(record_user_inputs_unit,'(E16.7,A)')TD,' quasi-peak discharge time constant, TD'
    write(post_process_info_unit,*)'quasi-peak discharge time constant, TD=',TD
  
    write(*,*)'Enter quasi-peak charging time constant, TC:'
    read(*,*)TC
    write(record_user_inputs_unit,'(E16.7,A)')TC,' quasi-peak charging time constant, TC'
    write(post_process_info_unit,*)'quasi-peak charging time constant, TC=',TC
  
  end if

! The frequency step is related to the filter resolution bandwidth

  FFILTER_fmin=-2d0*RBW
  FFILTER_fmax=2d0*RBW
  
! get the period of time domain data
  tmax=dt*function_of_time(time_function_number)%n_timesteps
    
  write(*,*)'tmax=',tmax,' s'
  write(post_process_info_unit,*)'tmax=',tmax,' s'
  
  FFILTER_fstep=1d0/tmax
  
  FFILTER_nf=FFILTER_fmax/FFILTER_fstep                ! number of frequency samples in 2*RBW

  write(*,*)
  write(*,*)'Set up the band pass filter:'

  FILTER_c=RBW/(2d0*sqrt(log(2d0)))                  ! Gaussian width

  write(*,*)'EMI filter RBW  =',real(RBW),' Hz'
  write(*,*)'EMI filter fmin =',real(FFILTER_fmin),' Hz'
  write(*,*)'EMI filter fmax =',real(FFILTER_fmax),' Hz'
  write(*,*)'EMI filter fstep=',real(FFILTER_fstep),' Hz'
  write(*,*)'EMI filter nf   =',FFILTER_nf
  write(*,*)'FILTER_c       =',real(FILTER_c),' Hz'
  write(*,*)

  write(post_process_info_unit,*)
  write(post_process_info_unit,*)'EMI filter RBW  =',real(RBW),' Hz'
  write(post_process_info_unit,*)'EMI filter fmin =',real(FFILTER_fmin),' Hz'
  write(post_process_info_unit,*)'EMI filter fmax =',real(FFILTER_fmax),' Hz'
  write(post_process_info_unit,*)'EMI filter fstep=',real(FFILTER_fstep),' Hz'
  write(post_process_info_unit,*)'EMI filter nf   =',FFILTER_nf
  write(post_process_info_unit,*)'FILTER_c       =',real(FILTER_c),' Hz'
  write(post_process_info_unit,*)

  ALLOCATE( filter_response(-FFILTER_nf:FFILTER_nf) )
  ALLOCATE( filtered_signal(-FFILTER_nf:FFILTER_nf) )
  ALLOCATE( Amag(-FFILTER_nf:FFILTER_nf) )
  ALLOCATE( Aphase(-FFILTER_nf:FFILTER_nf) )
  ALLOCATE( Aw(-FFILTER_nf:FFILTER_nf) )

  do floop=-FFILTER_nf,FFILTER_nf
    f=dble(floop)*FFILTER_fstep
    filter_response(floop)=exp(-(f/FILTER_c)**2 ) 
  end do

! Set the Fourier transform parameters

  fmin=FSWEEP_fmin+FFILTER_fmin   ! sufficient range for frequency sweep and filter bandwidth
  fmax=FSWEEP_fmax+FFILTER_fmax   ! sufficient range for frequency sweep and filter bandwidth
  fstep=FFILTER_fstep
  n_frequencies=1+NINT((fmax-fmin)/fstep)
  
  write(*,*)'Fourier transform parameters:'
  write(*,*)'fmin  =',fmin,' Hz'
  write(*,*)'fmax  =',fmax,' Hz'
  write(*,*)'fstep =',fstep,' Hz'
  write(*,*)'n_frequencies==',n_frequencies
  write(post_process_info_unit,*)'Fourier transform parameters:'
  write(post_process_info_unit,*)'fmin  =',fmin,' Hz'
  write(post_process_info_unit,*)'fmax  =',fmax,' Hz'
  write(post_process_info_unit,*)'fstep =',fstep,' Hz'
  write(post_process_info_unit,*)'n_frequencies=',n_frequencies

! Set up the envelope function parameters

  write(*,*)
  write(*,*)'Set up the envelope function parameters:'

! The filtered signal is periodic with period tmax_envelope
  tmax_envelope=1d0/FFILTER_fstep

! The maximum frequency of the envelope function is FFILTER_fmax.
! Choose n_envelope so we have say 20 samples in each period of this max frequency

  n_envelope=NINT(20d0*tmax_envelope*FFILTER_fmax)

  dt_envelope=tmax_envelope/dble(n_envelope-1)

  write(*,*)
  write(*,*)'Tmax_envelope=',real(tmax_envelope),' s'
  write(*,*)'dt_envelope  =',real(dt_envelope),' s'
  write(*,*)'n_envelope   =',n_envelope
  write(post_process_info_unit,*)
  write(post_process_info_unit,*)'Tmax_envelope=',real(tmax_envelope),' s'
  write(post_process_info_unit,*)'dt_envelope  =',real(dt_envelope),' s'
  write(post_process_info_unit,*)'n_envelope   =',n_envelope

  ALLOCATE( envelope(n_envelope) )

! Generate the frequency list
  function_of_frequency(freq_function_number)%n_frequencies=n_frequencies
  
  ALLOCATE ( function_of_frequency(freq_function_number)%frequency(1:n_frequencies) ) 
  do frequency_loop=1,n_frequencies     
    function_of_frequency(freq_function_number)%frequency(frequency_loop)=fmin+(frequency_loop-1)*fstep 
  end do

! Allocate Fourier transform data  
  ALLOCATE ( function_of_frequency(freq_function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(freq_function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(freq_function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(freq_function_number)%dB(1:n_frequencies) )

! calculate the Fourier integral at each of the frequencies specified
  do frequency_loop=1,n_frequencies
    
    w=2d0*pi*function_of_frequency(freq_function_number)%frequency(frequency_loop)
    
    integral=(0d0,0d0)

    do timestep=1,function_of_time(time_function_number)%n_timesteps
    
      integral=integral+function_of_time(time_function_number)%value(timestep)	&
                       *exp(-j*w*function_of_time(time_function_number)%time(timestep))*dt
		       
    end do ! next timestep

! Set the value at this frequency with scaling appropriate for RMS output for EMI receiver

    function_of_frequency(freq_function_number)%value(frequency_loop)=integral*sqrt(2d0)/tmax
    
    function_of_frequency(freq_function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(freq_function_number)%value(frequency_loop))
    function_of_frequency(freq_function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(freq_function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(freq_function_number)%value(frequency_loop))   )
    function_of_frequency(freq_function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(freq_function_number)%magnitude(frequency_loop))

  end do ! next frequency 
  
! Write the raw Fourier transform data to file
  
  CALL write_Frequency_Domain_Data(freq_function_number)

  ! Allocate frequency domain EMI response data

! Generate the frequency list for the EMI reciever sweep
  function_of_frequency(EMI_function_number)%n_frequencies=n_FSWEEP
  
  ALLOCATE ( function_of_frequency(EMI_function_number)%frequency(1:n_FSWEEP) )
  do frequency_loop=1,n_FSWEEP      
    function_of_frequency(EMI_function_number)%frequency(frequency_loop)=FSWEEP_fmin+(frequency_loop-1)*FSWEEP_fstep   
  end do
  ALLOCATE ( function_of_frequency(EMI_function_number)%value(1:n_FSWEEP) )
  ALLOCATE ( function_of_frequency(EMI_function_number)%magnitude(1:n_FSWEEP) )
  ALLOCATE ( function_of_frequency(EMI_function_number)%phase(1:n_FSWEEP) )
  ALLOCATE ( function_of_frequency(EMI_function_number)%dB(1:n_FSWEEP) )
 
  do frequency_loop=1,n_FSWEEP

    f=function_of_frequency(EMI_function_number)%frequency(frequency_loop)    ! Central frequency of filter
  
! get the frequency sample corresponding to this frequency

    nf0=NINT( 1d0+(f-fmin)/fstep )

! Apply the IF filter function

! get the minimum and maximum frequency samples in the IF filter band
  
    nfm=nf0-FFILTER_nf
    nfp=nf0+FFILTER_nf
    
    if (nfm.LT.1) then
      write(*,*)'ERROR: Filter_fmin<0'
      STOP 1
    end if
  
    if (nfp.GT.n_frequencies) then
      write(*,*)'ERROR: Filter_fmax>FFT_fmax'
      STOP 1
    end if

!  Multiply the siganl frequency response by the filter frequency response 
    Amag(:)=0d0
    Aphase(:)=0d0
    Aw(:)=0d0
      
    do i=nfm,nfp             ! i is the FFT sample number
  
      nfsample=i-nf0         ! nfsample is the filter sample number

      filtered_signal(nfsample)=filter_response(nfsample)*function_of_frequency(freq_function_number)%value(i)
    
      Amag(nfsample)=abs(filtered_signal(nfsample))
      Aphase(nfsample)=atan2(imag(filtered_signal(nfsample)),dble(filtered_signal(nfsample)))
      Aw(nfsample)=2d0*pi*(fmin+(i-1)*fstep)
    
    end do

! Calculate the envelope function

    Aw0=Aw(0)               ! central frequency of the filter
  
    do tloop=1,n_envelope
  
      t=dble(tloop-1)*dt_envelope

! CJS form using equation (7) which has a single sum and should be more efficient

      envelope(tloop)=0d0
    
      CT=0d0
      ST=0d0
      do ii=-FFILTER_nf,FFILTER_nf
        CT=CT+Amag(ii)*cos((Aw(ii)-Aw0)*t + Aphase(ii))
        ST=ST+Amag(ii)*sin((Aw(ii)-Aw0)*t + Aphase(ii))
      end do  
    
      envelope(tloop)=envelope(tloop)+sqrt(CT*CT+ST*ST)
     
    end do ! next timestep
    
    Vpeak=0d0
    do tloop=1,n_envelope
      Vpeak=max( Vpeak,envelope(tloop) )   
    end do  

    Vrms=0d0
    do tloop=1,n_envelope
      Vrms=Vrms+abs(envelope(tloop))**2    
    end do
    Vrms=sqrt(Vrms/dble(n_envelope))

    Vavg=0d0
    do tloop=1,n_envelope
      Vavg=Vavg+abs(envelope(tloop))   
    end do
    Vavg=Vavg/dble(n_envelope)
           
    if (detector_type.EQ.detector_type_peak) then
      Vout=Vpeak      
    elseif (detector_type.EQ.detector_type_rms) then
      Vout=Vrms      
    elseif (detector_type.EQ.detector_type_avg) then
      Vout=Vavg      
    elseif (detector_type.EQ.detector_type_quasi_peak) then
  
! Crude search for Vqp. Could be done much more efficiently I think...

! initial range of search
      Vqpmin=Vavg
      Vqpmax=Vpeak
    
      do i=1,10    ! number of search loops

        Vqp=(Vqpmax+Vqpmin)/2d0

! Calculate the RHS
        RHS=Vqp*n_envelope/TD
  
 ! Calculate the LHS
        LHS=0d0
        do tloop=1,n_envelope
          if (envelope(tloop).GT.Vqp) then
            LHS=LHS+abs(envelope(tloop)-Vqp)/TC   
          end if
        end do
     
        if (LHS.GT.RHS) then
! The trial value of Vqp is too low so increase the minimum value of the range
          Vqpmin=Vqp
        else
! The trial value of Vqp is too high so decrease the maximum value of the range
          Vqpmax=Vqp
        end if
    
      end do
 
      Vquasi_peak=(Vqpmax+Vqpmin)/2d0
      Vout=Vquasi_peak
      
    else
    
      write(*,*)'Unknown detector_type:',detector_type
      STOP 1
    
    end if
    
    write(*,*)'frequency=',real(f),'Vpeak=',real(20.0*log10(1e6*Vpeak)), &
                                   ' Vrms=',real(20.0*log10(1e6*Vrms)),  &
                                   ' Vquasi-peak=',real(20.0*log10(1e6*Vquasi_peak)), &
                                   ' Vavg=',real(20.0*log10(1e6*Vavg))

    function_of_frequency(EMI_function_number)%value(frequency_loop)=Vout
    
    function_of_frequency(EMI_function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(EMI_function_number)%value(frequency_loop))
    function_of_frequency(EMI_function_number)%phase(frequency_loop)=0d0
    function_of_frequency(EMI_function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(EMI_function_number)%magnitude(frequency_loop))
 
  end do ! next frequency in the specified frequency sweep

! Write the EMI receiver data to file
  
  CALL write_Frequency_Domain_Data(EMI_function_number)
  

  DEALLOCATE( envelope )
  DEALLOCATE( filter_response )
  DEALLOCATE( filtered_signal )
  DEALLOCATE( Amag )
  DEALLOCATE( Aphase )
  DEALLOCATE( Aw )

  CALL Deallocate_post_data()

  RETURN
  
  
END SUBROUTINE EMI_receiver_model
