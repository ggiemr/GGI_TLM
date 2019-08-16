PROGRAM create_test_datasets

IMPLICIT NONE

integer :: ncycles,ns,nh
real*8    :: f0,fh,period
integer :: nharmonics
real*8    :: a50,a0,an

character(LEN=256) :: filename

integer :: tot_nsamples

real*8    :: time,dt
real*8    :: a
integer :: i,ih,ic

real*8,allocatable :: V(:)

real*8 :: pi

integer*4 idnum

! function types

real*8 urand
real*8 nrand

! START

write(*,*)'Enter the square wave frequency, f0'
read(*,*)f0

write(*,*)'Enter the square wave amplitude'
read(*,*)a0

write(*,*)'Enter the number of samples in one cycle of the square wave'
read(*,*)ns

write(*,*)'Enter the number of cycles of the square wave'
read(*,*)ncycles

write(*,*)'Enter the number of odd harmonics'
read(*,*)nh

write(*,*)'Enter the amplitude of 50Hz signal to add'
read(*,*)a50

write(*,*)'Enter the amplitude of gaussian noise to add'
read(*,*)an

write(*,*)'Enter the filename for the signal output'
read(*,'(A)')filename
  
idnum=1347943

pi=4d0*atan(1d0)
period=1d0/f0
dt=period/ns

tot_nsamples=ns*ncycles

ALLOCATE( V(1:ns+1) )

! Build the waveform for one cycle

time=0d0
V(:)=0d0

do i=1,ns+1

  time=(i-1)*dt

  do ih=1,nh
  
    a=a0*4d0/(pi*(2*ih-1))
    
    V(i)=V(i)+a*sin((2*ih-1)*2d0*pi*f0*time)
      
  end do ! next odd harmonic

end do ! next sample

open(unit=10,file=trim(filename))

do ic=1,ncycles

  do i=1,ns

    time=((ic-1)*ns+(i-1))*dt

    write(10,*)time,V(i)+a50*sin(2d0*pi*50d0*time)+an*nrand(idnum)
  
  end do ! next sample
  
end do ! next cycle

! write a final timestep
time=time+dt
write(10,*)time,V(ns+1)+a50*sin(2d0*pi*50d0*time)

close(unit=10)

DEALLOCATE( V )

open(unit=12,file='square_wave_fourier_series.fout')

do ih=1,nh

  a=a0*4d0/(pi*(2*ih-1))
  
! write frequency, Vrms, VdBm, VdBuV
  write(12,1000)(2*ih-1)*f0,a/sqrt(2d0),10d0*log10((a*a/(2d0*50d0))/1000d0),20d0*log10(1d6*a/(sqrt(2d0)))
1000 format(4ES16.6)

end do ! next odd harmonic

close(unit=12)

open(unit=14,file='input_signal_summary.dat')

write(14,*)'Number of samples =',tot_nsamples
write(14,'(A,ES16.6)')'Timestep =',dt
write(14,'(A,ES16.6)')'F_Nyquist =',1d0/dt
write(14,'(A,ES16.6)')'F_max =',1d0/(2d0*dt)
write(14,'(A,ES16.6)')'Min bandwidth =',1d0/(tot_nsamples*dt)

close(unit=14)

END PROGRAM create_test_datasets
!
! NAME
!     urand
!
! DESCRIPTION
!     return a random real from a uniform distribution
!     Code from Numerical Recipes
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 23/09/2013 CJS
!
!
       
       function urand(idnum)

IMPLICIT NONE
       
       real*8 urand
       integer*4 idnum
       
       integer*4 ia,im,iq,ir
       real*8 am
       real*8 k
       
       ia=16807
       im=2147483647
       am=1.0/im
       iq=127773
       ir=2836
       
       k=idnum/iq
       idnum=ia*(idnum-k*iq)-ir*k
       if (idnum.lt.0) idnum=idnum+im
       urand=am*idnum
       
       return
       end
!
! NAME
!     nrand
!
! DESCRIPTION
!     return a random real from a normal distribution
!     Code from Numerical Recipes
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 23/09/2013 CJS
!
!
       function nrand(idnum)

IMPLICIT NONE
       
       real*8 nrand
       integer*4 idnum
       
       real*8 v1,v2,rsqr,fac
       
       real*8 urand
       
10     continue       
         v1=2d0*urand(idnum)-1d0
         v2=2d0*urand(idnum)-1d0
         rsqr=v1*v1+v2*v2
       if ( (rsqr.ge.1d0).OR.(rsqr.eq.0d0) ) goto 10
       fac=sqrt(-2d0*log(rsqr)/rsqr)
       nrand=v2*fac
              
       return
       end
