
PROGRAM wire_over_lossy_ground

IMPLICIT NONE

real*8 :: fmin,fmax,f,w
real*8 :: logfmin,logfmax,logfstep,logf

real*8 :: sigma,a,h,epsr

real*8 :: test1,test2,test3,test4,test5

complex*16 :: Z,Y,gamma,beta
complex*16 :: Zc,Yc,gammac,betac
complex*16 :: Zc2,Yc2,gammac2,betac2
complex*16 :: Zw,Yw,gammaw,betaw
complex*16 :: delta,k1,k2,ka,Jc,ALPHA,QmjP,NmjM

integer :: nf,floop

! complex plane search stuff for Wait's solution

integer    :: itest,min_test_point
real*8     :: dbeta_r,dbeta_i
complex*16 :: beta_test(0:4),new_beta,beta_save
real*8     :: error_test(0:4),error,min_error,initial_error,final_error,spatial_avg_final_error

real*8     :: prop_dist

integer    :: iteration
integer    :: max_iterations=1000

! Wait model mode search stuff

integer,parameter :: max_modes=10000

integer :: n_trial_modes,n_modes,mode

complex*16 :: beta_err_min(max_modes)
complex*16 :: beta_found(max_modes)

integer    :: n_re_beta,ire
real*8     :: re_beta_min,re_beta_max,d_re_beta,re_beta

integer    :: n_im_beta,iim
real*8     :: im_beta_min,im_beta_max,d_im_beta,im_beta

real*8     :: beta_diff,min_beta_diff
integer    :: mode_test
logical    :: already_found_this_mode

logical :: plot_beta_map

real*8,allocatable     :: err_array(:,:)
complex*16,allocatable :: beta_array(:,:)

character(LEN=80) :: FH2_DIRECTORY

! constants

real*8,parameter :: small=1d-12

real*8,parameter :: large=1d12

real*8,parameter :: pi=3.141592653589793d0

real*8,parameter :: Z0=376.73031346177d0

real*8,parameter :: Y0=1D0/Z0

complex*16,parameter :: j=(0d0,1d0)

real*8,parameter :: eps0=8.8541878176D-12

real*8,parameter :: mu0=3.141592653589793d0*4D-7

real*8,parameter :: c0=2.99792458D8


! START

! Set up default test case values

a=0.013      ! wire radius
h=20.0       ! wire height over ground
sigma=0.001  ! ground electrical conductivity
epsr=1.0     ! ground relative permittivity

fmin=1d4
fmax=10d6
nf=101

prop_dist=1000d0

n_re_beta=101
n_im_beta=101

write(*,*)'GGI_wire_over_lossy_ground'
write(*,*)'Calculate the complex propagation constant of modes on a wire over lossy ground'
write(*,*)'using the following solution methods:'
write(*,*)'Full method due to Carson'
write(*,*)'Approximate method due to Carson'
write(*,*)'Rigorous method due to Wait (multi-mode)'
write(*,*)'This software uses the Fortran complex bessel function module mod_zbes.f90'
write(*,*)'downloaded from: dl.acm.org/doi/10.1145/1916461.1916471'
write(*,*)
write(*,*)'Enter the wire radius (m)'
read(*,*)a
write(*,*)'Enter the wire height over the ground (m)'
read(*,*)h
write(*,*)'Enter the ground electrical conductivity (S/m)'
read(*,*)sigma
write(*,*)'Enter the ground relative permittivity'
read(*,*)epsr
write(*,*)'Enter the minimum frequency for analysis (Hz)'
read(*,*)fmin
write(*,*)'Enter the maximum frequency for analysis (Hz)'
read(*,*)fmax
write(*,*)'Enter the number of frequencies to evaluate (log frequency scale)'
read(*,*)nf

write(*,*)'Writing configuration to file: wire_over_lossy_ground_configuration.dat'

open(unit=10,file='wire_over_lossy_ground_configuration.dat')

write(10,*)real(a),'  wire radius (m)'
write(10,*)real(h),'  wire height over the ground (m)'
write(10,*)real(sigma),'  ground electrical conductivity (S/m)'
write(10,*)real(epsr),'  ground relative permittivity'
write(10,*)real(fmin),'  minimum frequency for analysis'
write(10,*)real(fmax),'  maximum frequency for analysis'
write(10,*)nf,'  number of frequencies to evaluate (log frequency scale)'

close(unit=10)

ALLOCATE( err_array(1:n_re_beta,1:n_im_beta) )
ALLOCATE( beta_array(1:n_re_beta,1:n_im_beta) )

open(unit=10,file='k0.out')
open(unit=12,file='Carson_wire_over_lossy_ground.out')
open(unit=14,file='Approx_Carson_wire_over_lossy_ground.out')
open(unit=16,file='Wait_wire_over_lossy_ground.out')
open(unit=18,file='Wait_wire_over_lossy_ground_all_modes.out')

logfmin=log10(fmin)
logfmax=log10(fmax)
if (nf.NE.1) then
  logfstep=(logfmax-logfmin)/dble(nf-1)
  plot_beta_map=.FALSE.
else
  logfstep=0d0
  plot_beta_map=.TRUE.
end if

do floop=1,nf

  logf=logfmin+dble(floop-1)*logfstep
  f=10.0**(logf)
  w=2d0*pi*f
  
  write(*,*)
  write(*,*)'Frequency=',real(f)
  
  write(10,'(2ES12.4)')f,w*sqrt(mu0*eps0)
    
! Approximate Carson model (no magnetic materials)

  delta=sqrt(2.0/(w*mu0*sigma))
  Zc=w*mu0/8.0+(j*w*mu0/(2.0*pi))*( log(sqrt(2.0)*delta/a+log(0.926) ) )
  Yc=2.0*pi*j*w*eps0/log(2.0*h/a)    ! admittance for wire over perfect ground plane
  gammac=sqrt(Zc*Yc)
  k2=w*sqrt(mu0*eps0*epsr)
  
  write(*,*)'Zc =',cmplx(Zc),' Yc =',cmplx(Yc),' gammac =',cmplx(gammac)
    
  test1=real(2.0*k2*h)
  betac=gammac/j
  write(12,'(4ES12.4)')f,real(gammac),imag(gammac),20*log10(abs(exp(-gammac*prop_dist)))
  
  
! Second Approximate Carson model 

  k2=w*sqrt(mu0*(eps0*epsr-j*sigma/w))
  CALL calc_carson_integral(k2,f,a,h,sigma,epsr,Jc)

  Zc2=(j*w*mu0/(2.0*pi))*(log(2.0*h/a)+Jc)
  Yc2=2.0*pi*j*w*eps0/log(2.0*h/a)
  k1=w*sqrt(mu0*eps0)
  betac2=k1*sqrt(1d0+Jc/log(2.0*h/a))
  
  gammac2=j*betac2
  write(*,*)'Zc2=',cmplx(Zc2),' Yc2=',cmplx(Yc2),' gammac2=',cmplx(gammac2)
  gammac2=sqrt(Zc2*Yc2)  
  write(*,*)'Zc2=',cmplx(Zc2),' Yc2=',cmplx(Yc2),' gammac2=',cmplx(gammac2)
      
  write(14,'(4ES12.4)')f,real(gammac2),imag(gammac2),20*log10(abs(exp(-gammac2*prop_dist)))
    
! Wait model: starting beta from the Carson model 

  ka=k1*k2/sqrt(k2*k2+k1*k1)
  if(imag(ka).GT.0d0) ka=-ka

  write(*,*)
  write(*,*)'Carson  model      , beta=',betac
  write(*,*)'Approx Carson model, beta=',betac2
  write(*,*)
  write(*,*)'Branch points:'
  write(*,*)'beta=k0=',k1
  write(*,*)'beta=k2=',k2
  write(*,*)'beta=ka=',ka
  write(*,*)"Calculating resuts of Wait's model:"
  
  err_array(:,:)=0d0
  beta_array(:,:)=(0d0,0d0)
  
! Set search range for beta

  re_beta_min=k1*0.95
  re_beta_max=k1*1.1
  d_re_beta=(re_beta_max-re_beta_min)/dble(n_re_beta-1)
  
  write(*,*)'re_beta_min=',re_beta_min,' re_beta_max=',re_beta_max,' d_re_beta',d_re_beta

  im_beta_min=imag(sqrt(Zc2*Yc2)/j) *5d0                  ! range derived from Carson model
  im_beta_min=min(im_beta_min,-0.5e-3)
  
  im_beta_max=0d0
  d_im_beta=(im_beta_max-im_beta_min)/dble(n_im_beta-1)
  
  write(*,*)'im_beta_min=',im_beta_min,' im_beta_max=',im_beta_max,' d_im_beta',d_im_beta
  
  min_beta_diff=sqrt(d_re_beta**2+d_im_beta**2)/10.0   ! measure for testing if modes are equivalent
  
! evaluate error on grid

  write(*,*)'evaluate error on grid'
  
  if (plot_beta_map) then
    open(unit=50,file='beta_err_map.dat')
  end if
  
  do ire=1,n_re_beta
  
    re_beta=re_beta_min+dble(ire-1)*d_re_beta

    do iim=1,n_im_beta
  
      im_beta=im_beta_min+dble(iim-1)*d_im_beta
      
      beta=dcmplx(re_beta,im_beta)
            
      beta_array(ire,iim)=beta
      
      INCLUDE 'calc_error_Wait.F90'  
        
      err_array(ire,iim)=error
      
      if (plot_beta_map) then
        write(50,'(3ES12.4)')real(beta_array(ire,iim)),imag(beta_array(ire,iim)),real(log10(err_array(ire,iim)))
      end if
      
    end do
    
    if (plot_beta_map) then
      write(50,*)
    end if
    
  end do
  
  if (plot_beta_map) then
    close(unit=50)
  end if

! find local minima in error

  n_trial_modes=0
  
  do iim=2,n_im_beta-1
  
    do ire=2,n_re_beta-1
    
      error=err_array(ire,iim)
          
      if ( (error.LT.err_array(ire+1,iim  )).AND.   &
           (error.LT.err_array(ire-1,iim  )).AND.   &
           (error.LT.err_array(ire  ,iim+1)).AND.   &
           (error.LT.err_array(ire  ,iim-1))      ) then

! found a local minimum indicating a mode

        n_trial_modes=n_trial_modes+1 
        beta_err_min(n_trial_modes)=beta_array(ire,iim)
           
      end if
      
    end do
    
  end do

  write(16,'(ES12.4)',ADVANCE='NO')f

! loop over local minima in error and accurately calculate the mode

  n_modes=0
  
  do mode=1,n_trial_modes

! set starting point for mode search

    beta=beta_err_min(mode)
  
    k1=w*sqrt(mu0*eps0)
    k2=w*sqrt(mu0*(eps0*epsr-j*sigma/w))
  
    beta_test(0)=beta
    
    iteration=0
    initial_error=0d0
  
    dbeta_r=dble(beta)/100d0
    dbeta_i=dble(beta)/100d0
  
5000 CONTINUE   ! start of a new loop of the optimisation process

    iteration=iteration+1

! set the points of the simplex
    beta_test(1)=beta_test(0)-dbeta_r
    beta_test(2)=beta_test(0)+dbeta_r
    beta_test(3)=beta_test(0)-j*dbeta_i
    beta_test(4)=beta_test(0)+j*dbeta_i
  
    min_error=1e30
    min_test_point=-1
  
    do itest=0,4
  
      beta=beta_test(itest)
      INCLUDE 'calc_error_Wait.F90'    
      error_test(itest)=error
          
      if (iteration.EQ.1) then                   ! NEW: include all 5 points in initial error calculation
        initial_error=initial_error+error/5d0
      end if
    
      if (error.LT.min_error) then
        min_error=error
        min_test_point=itest
      end if
  
    end do
  
    if (min_test_point.EQ.-1) then
      write(*,*)'ERROR evaulating error in Wait solution:'
      STOP 1
    end if
  
    if (min_test_point.EQ.0) then
!   reduce the size of the simplex
      dbeta_r=dbeta_r/2d0
      dbeta_i=dbeta_i/2d0
    else
!   put the minimum error solution at the centre of the simplex  
      beta_test(0)=beta_test(min_test_point)
    end if
  
    if (iteration.LT.max_iterations) GOTO 5000
 
! recreate the simplex around the minimum point 
    betaw=beta_test(0)
    
! set the points of the simplex surrounding the minimum point but increase the simplex size to 
! detect minima near branch lines in the solution which are not true solutions

    dbeta_r=dble(beta)/1000d0
    dbeta_i=dble(beta)/1000d0
    
    beta_test(1)=beta_test(0)-dbeta_r
    beta_test(2)=beta_test(0)+dbeta_r
    beta_test(3)=beta_test(0)-j*dbeta_i
    beta_test(4)=beta_test(0)+j*dbeta_i

! evaluate the error on all points of the simplex and calculate the final error

    spatial_avg_final_error=0d0
    do itest=0,4  
      beta=beta_test(itest)
      INCLUDE 'calc_error_Wait.F90'    
      error_test(itest)=error
      spatial_avg_final_error=spatial_avg_final_error+error/5d0
    end do
    
    final_error=error_test(0)
       
    if ( (real(final_error/initial_error).LT.1d-6).AND. &
         (real(spatial_avg_final_error/initial_error).LT.1d-1).AND. &
         (real(betaw).GE.0d0).AND. (imag(betaw).LT.0d0)  ) then
         
      write(*,'(A,I4,A,F10.6,A,F10.6,A,ES10.2,A,ES10.2)')'    mode=',mode, &
              ' betaw =',real(betaw),'+j',imag(betaw),'  point_err=',real(final_error/initial_error), &
              ' avg_err=',real(spatial_avg_final_error/initial_error)
         
! assume that the error is tending to zero and this is a real forward propagating mode with attenuation
    
      if (mode.EQ.1) then  ! this must be a new mode
        n_modes=1
        beta_found(n_modes)=betaw
        write(16,'(3ES12.4)',ADVANCE='NO')real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
        write(18,'(4ES12.4)')f,real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
         
        write(*,*)
        write(*,*)'mode=',n_modes
        write(*,*)'Z=',Zw
        write(*,*)'Y=',Yw
        write(*,*)'R=',real(Zw)
        write(*,*)'L=',real(Zw/(j*w))
        write(*,*)'G=',real(Yw)
        write(*,*)'C=',real(Yw/(j*w))
          
       else
    
! check to see if we have already found this mode
        already_found_this_mode=.FALSE.
        do mode_test=1,n_modes
          beta_diff=abs(betaw-beta_found(mode_test))
          if (beta_diff.LT.min_beta_diff) then
            already_found_this_mode=.TRUE.
          end if
        end do
    
        if(.NOT.already_found_this_mode) then
          n_modes=n_modes+1
          beta_found(n_modes)=betaw    
          write(16,'(3ES12.4)',ADVANCE='NO')real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
          write(18,'(4ES12.4)')f,real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
          
          
          write(*,*)
          write(*,*)'mode=',n_modes
          write(*,*)'Z=',Zw
          write(*,*)'Y=',Yw
          write(*,*)'R=',real(Zw)
          write(*,*)'L=',real(Zw/(j*w))
          write(*,*)'G=',real(Yw)
          write(*,*)'C=',real(Yw/(j*w))
                    
        end if
      
      end if  ! mode=1 or not...
      
    else
    
      write(*,'(A,I4,A,F10.6,A,F10.6,A,ES10.2,A,ES10.2)')'*** mode=',mode, &
              ' betaw =',real(betaw),'+j',imag(betaw),'  point_err=',real(final_error/initial_error), &
              ' avg_err=',real(spatial_avg_final_error/initial_error)

    end if    ! relative error is small
    
  end do ! next mode
  
  write(16,*)

8888 CONTINUE

  flush(unit=10)
  flush(unit=12)
  flush(unit=14)
  flush(unit=16)
  flush(unit=18)

end do

close(unit=10) 
close(unit=12) 
close(unit=14) 
close(unit=16)
close(unit=18)
 
DEALLOCATE( err_array )
DEALLOCATE( beta_array )

write(*,*)
write(*,*)'Output files:'
write(*,*)"file='k0.out' : Free space wave number"
write(*,*)"file='Carson_wire_over_lossy_ground.out' : Carson model"
write(*,*)"file='Approx_Carson_wire_over_lossy_ground.out' : Approximate Carson model"
write(*,*)"file='Wait_wire_over_lossy_ground.out' : Wait model, first two distinct modes"
write(*,*)"file='Wait_wire_over_lossy_ground_all_modes.out' : Wait model, all modes found"
write(*,*)
write(*,*)'Output file formats for wire_over_lossy_ground: single mode.'
write(*,*)' f(Hz) attenuation constant, propagation constant, attenuation over 1km'

write(*,*)
write(*,*)'Re-run with this configuration with the command: '
write(*,*)'wire_over_lossy_ground < wire_over_lossy_ground_configuration.dat'

STOP

END PROGRAM wire_over_lossy_ground
!
! ________________________________________________________________________________
!
!
SUBROUTINE calc_carson_integral(k2,f,a,h,sigma,epsr,Jc)

IMPLICIT NONE

complex*16 :: k2
real*8     :: f,a,h,sigma,epsr
complex*16 :: Jc

! local variables

real*8 :: l,lmax,dl
real*8 :: w
real*8 :: fmin
complex*16 :: u,gammag
complex*16 :: integrand
integer :: nl,i

character(LEN=80) :: command
integer :: exit_stat

real*8,parameter :: pi=3.141592653589793d0

complex*16,parameter :: j=(0d0,1d0)

real*8,parameter :: eps0=8.8541878176D-12

real*8,parameter :: mu0=3.141592653589793d0*4D-7

! START

nl=1000
fmin=1d-4

!nl=10000
!fmin=1d-5

lmax=-log(fmin)/(2d0*h)
dl=lmax/dble(nl)

! From Olsen paper

Jc=(0d0,0d0)

do i=1,nl

  l=dble(i)*dl-dl/2d0
  u=sqrt(l*l-k2*k2)            ! sign of sqare root doesn't appear to make a difference
  integrand=-(2d0/(k2*k2))*((u-l)*exp(-2d0*l*h))   ! note sign chnge from Olsen paper to correct sign discrepancy
  Jc=Jc+dl*integrand
  
end do

RETURN

! The following code uses different forms for the equation for testing.

write(*,*)'Carson Jc integral (Olsen)  , f=',f,' Jc=',Jc

! from Papadopoulos paper

w=2d0*3.1415926535d0*f
gammag=sqrt(j*w*mu0*sigma)

Jc=(0d0,0d0)

do i=1,nl

  l=dble(i)*dl-dl/2d0
! Papadopoulos equation 4
  integrand=2d0*exp(-2d0*l*h)*cos(l*a)/(l+sqrt(l*l+gammag*gammag))
  Jc=Jc+dl*integrand
  
!  write(20,'(3ES12.4)')real(l),real(integrand),imag(integrand)

end do

write(*,*)'Carson Jc integral P1       , f=',f,' Jc=',Jc

gammag=sqrt(j*w*mu0*(sigma+j*w*eps0*epsr))

Jc=(0d0,0d0)

do i=1,nl

  l=dble(i)*dl-dl/2d0
! Papadopoulos equation 4
  integrand=2d0*exp(-2d0*l*h)*cos(l*a)/(l+sqrt(l*l+gammag*gammag))
  Jc=Jc+dl*integrand
  
!  write(20,'(3ES12.4)')real(l),real(integrand),imag(integrand)

end do

write(*,*)'Carson Jc integral P2(gamma), f=',f,' Jc=',Jc

gammag=sqrt(j*w*mu0*(sigma+j*w*eps0*epsr))

Jc=(0d0,0d0)

do i=1,nl

  l=dble(i)*dl-dl/2d0
! Papadopoulos equation 4
  integrand=2d0*exp(-2d0*l*h)/(l+sqrt(l*l+gammag*gammag))
  Jc=Jc+dl*integrand
  
!  write(20,'(3ES12.4)')real(l),real(integrand),imag(integrand)

end do

write(*,*)'Carson Jc integral P3 (cos) , f=',f,' Jc=',Jc

RETURN

END SUBROUTINE calc_carson_integral
!
! ________________________________________________________________________________
!
!
SUBROUTINE wiat_integral(k1,k2,f,a,h,beta,ALPHA,QmjP,NmjM)

IMPLICIT NONE

complex*16 :: k1,k2
real*8     :: f,a,h
complex*16 :: beta,ALPHA,QmjP,NmjM

! local variables

real*8 :: l,lmax,dl
real*8 :: fmin
complex*16 :: K0a,K0h,arg
complex*16 :: u1,u2
complex*16 :: integrand,integrand2
integer :: nl,i

character(LEN=80) :: command
integer :: exit_stat,cmd_stat

logical :: plot_integrand=.FALSE.

! START

nl=1000
fmin=1d-4

!nl=2000
!fmin=0.5d-4

lmax=-log(fmin)/(2d0*h)
dl=lmax/dble(nl)

arg=(0d0,1d0)*a*sqrt(k1*k1-beta*beta) ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0a)

arg=(0d0,1d0)*2d0*h*sqrt(k1*k1-beta*beta) ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0h)

ALPHA=K0a-K0h

QmjP=(0d0,0d0)
NmjM=(0d0,0d0)

if (plot_integrand) then
  open(unit=20,file='integrand.dat')
end if

do i=1,nl

  l=dble(i)*dl-dl/2d0
  u1=sqrt(l*l+beta*beta-k1*k1)
  u2=sqrt(l*l+beta*beta-k2*k2)
  integrand=(exp(-2d0*u1*h)/(u1+u2))*cos(l*a)
  QmjP=QmjP+dl*integrand
  
  integrand2=k1*k1*(exp(-2d0*u1*h)/(k2*k2*u1+k1*k1*u2))*cos(l*a)
  NmjM=NmjM+dl*integrand2
  
  if (plot_integrand) then
    write(20,'(5ES12.4)')real(l),real(integrand),imag(integrand),real(integrand2),imag(integrand2)
  end if
  
end do

if (plot_integrand) then
  write(*,*)'Wait integrals, f=',real(f),' ALPHA=',cmplx(ALPHA),' QmjP=',cmplx(QmjP),' NmjM=',cmplx(NmjM)

  close(unit=20)
  command='gnuplot i2plot.plt'
  CALL execute_command_line(command,EXITSTAT=exit_stat,CMDSTAT=cmd_stat)

  STOP 1
  
end if

END SUBROUTINE wiat_integral
!
! ________________________________________________________________________________
!
!
SUBROUTINE calc_modified_Bessel_function(znu,arg,K0)

USE mod_zbes

IMPLICIT NONE

complex*16 :: znu,arg,K0

! local variables

complex*16 :: zz,H1n,H2n

integer    :: info

! Evaluate modified Bessel function

  if (imag(arg).GE.0d0) then
    zz=arg*exp((0d0,-3.1415926535d0)/2d0)
    CALL hankel2(znu,zz,H2n,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znu
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0=((0d0,-3.1415926535d0)/2d0)*exp(znu*(0d0,-3.1415926535d0)/2d0)*H2n
  else
    zz=arg*exp((0d0,3.1415926535d0)/2d0)
    CALL hankel1(znu,zz,H1n,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znu
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0=((0d0,3.1415926535d0)/2d0)*exp(znu*(0d0,3.1415926535d0)/2d0)*H1n
  end if
  
  RETURN

END SUBROUTINE calc_modified_Bessel_function

