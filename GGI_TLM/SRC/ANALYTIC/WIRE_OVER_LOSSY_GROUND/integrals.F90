!
! ________________________________________________________________________________
!
!
SUBROUTINE calc_carson_integral(k2,f,a,h,sigma,epsr,Jc,nlambda,lambda,dlambda)

IMPLICIT NONE

complex*16 :: k2
real*8     :: f,a,h,sigma,epsr
complex*16 :: Jc
integer :: nlambda
real*8 :: lambda(0:nlambda)
real*8 :: dlambda(0:nlambda)

! local variables

real*8 :: l
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

!open(unit=20,file='integrand.dat')

! From Olsen paper

Jc=(0d0,0d0)

do i=0,nlambda

  l=lambda(i)
  u=sqrt(l*l-k2*k2)            ! sign of sqare root doesn't appear to make a difference
  integrand=-(2d0/(k2*k2))*((u-l)*exp(-2d0*l*h))   ! note sign chnge from Olsen paper to correct sign discrepancy
  Jc=Jc+integrand*dlambda(i)
  
!  write(20,'(3ES12.4)')real(l),real(integrand),imag(integrand)

end do

write(*,*)'Carson Jc integral (Olsen)  , f=',f,' Jc=',Jc

RETURN


END SUBROUTINE calc_carson_integral
!
! ________________________________________________________________________________
!
!
SUBROUTINE Wait_integral(k1,k2,f,a,h,beta,ALPHA,QmjP,NmjM,plot_integrand,show_plots,Jc,nlambda,lambda,dlambda)

IMPLICIT NONE

complex*16 :: k1,k2
real*8     :: f,a,h
complex*16 :: beta,ALPHA,QmjP,NmjM,Jc
logical    :: plot_integrand,show_plots
integer :: nlambda
real*8 :: lambda(0:nlambda)
real*8 :: dlambda(0:nlambda)

! local variables

real*8 :: l
real*8 :: fmin
complex*16 :: K0a,K0h,kr,arg
complex*16 :: u1,u2
complex*16 :: integrand,integrand2

complex*16 :: exp_term

integer :: nl,i

character(LEN=80) :: command
integer :: exit_stat,cmd_stat

! START

kr=sqrt(k1*k1-beta*beta)
if(imag(kr).GT.0d0) kr=kr*(-1d0)   ! sign check for lossy propagation away from the wire

arg=(0d0,1d0)*a*kr ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0a)

!arg=(0d0,1d0)*2d0*h*kr ! argument of modified Bessel function
arg=(0d0,1d0)*sqrt(4d0*h*h+a*a)*kr  ! From original Wait paper
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0h)

ALPHA=K0a-K0h

QmjP=(0d0,0d0)
NmjM=(0d0,0d0)

if (plot_integrand) then
  open(unit=20,file='integrand.dat')
end if

do i=0,nlambda

  l=lambda(i)
  u1=sqrt(l*l+beta*beta-k1*k1)
  if(dble(u1).LT.0d0) u1=u1*(-1d0)   ! sign check. Olsen paper section III
  
  u2=sqrt(l*l+beta*beta-k2*k2)
  if(dble(u2).LT.0d0) u2=u2*(-1d0)   ! sign check. Olsen paper section III
  
  exp_term=exp(-u1*sqrt(4d0*h*h+a*a))     ! From original Wait paper
   
  integrand=(exp_term/(u1+u2))*cos(l*a)   ! From original Wait paper
  QmjP=QmjP+integrand*dlambda(i)
  
  integrand2=k1*k1*(exp_term/(k2*k2*u1+k1*k1*u2))*cos(l*a)   ! adapted From Olsen paper
  NmjM=NmjM+integrand2*dlambda(i)
  
  if (plot_integrand) then
    write(20,'(5ES12.4)')real(l),real(integrand),imag(integrand),real(integrand2),imag(integrand2)
  end if
  
end do

if (plot_integrand) then
  write(*,*)'Wait integrals, f=',real(f),' ALPHA=',cmplx(ALPHA)
  write(*,*)'Carson Approximation: ALPHA~',cmplx(log(2d0*h/a))
  write(*,*)' QmjP =',cmplx(QmjP),' NmjM=',cmplx(NmjM)
  write(*,*)'Carson~',cmplx(Jc/2d0),' NmjM~',cmplx(0d0)
  
  write(*,*)'sqrt(k1*k1-beta*beta)=',sqrt(k1*k1-beta*beta)
  
  close(unit=20)
  
  if(show_plots) then
    command='gnuplot i2plot.plt'
    CALL execute_command_line(command,EXITSTAT=exit_stat,CMDSTAT=cmd_stat)
  end if
  
!  STOP 1
  
end if

END SUBROUTINE Wait_integral
!
! ________________________________________________________________________________
!
!
SUBROUTINE Wait_integral2(k1,k2,f,a,h,beta,Kterm,Integral_term,ALPHA,QmjP,NmjM, &
                          plot_integrand,show_plots,Jc,nlambda,lambda,dlambda)

IMPLICIT NONE

complex*16 :: k1,k2
real*8     :: f,a,h
complex*16 :: beta,Kterm,Integral_term,ALPHA,QmjP,NmjM,Jc
logical    :: plot_integrand,show_plots
integer :: nlambda
real*8 :: lambda(0:nlambda)
real*8 :: dlambda(0:nlambda)

! local variables

real*8 :: l
real*8 :: fmin
complex*16 :: K0a,K0h,kr,arg
complex*16 :: u1,u2
complex*16 :: integrand,integrand2,integrand3

complex*16 :: exp_term

integer :: nl,i

character(LEN=80) :: command
integer :: exit_stat,cmd_stat

! START


kr=sqrt(k1*k1-beta*beta)
if(imag(kr).GT.0d0) kr=kr*(-1d0)   ! sign check for lossy propagation away from the wire

arg=(0d0,1d0)*a*kr ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0a)

!arg=(0d0,1d0)*2d0*h*kr ! argument of modified Bessel function
arg=(0d0,1d0)*sqrt(4d0*h*h+a*a)*kr  ! From original Wait paper
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0h)

ALPHA=K0a-K0h

Kterm=(1d0-beta*beta/(k1*k1))*(K0a-K0h)    ! Wait, equation 17

QmjP=(0d0,0d0)
NmjM=(0d0,0d0)
Integral_term=(0d0,0d0)

if (plot_integrand) then
  open(unit=20,file='integrand.dat')
end if

do i=0,nlambda

  l=lambda(i)
  u1=sqrt(l*l+beta*beta-k1*k1)
  if(dble(u1).LT.0d0) u1=u1*(-1d0)   ! sign check. Olsen paper section III
  
  u2=sqrt(l*l+beta*beta-k2*k2)
  if(dble(u2).LT.0d0) u2=u2*(-1d0)   ! sign check. Olsen paper section III
  
  exp_term=exp(-u1*sqrt(4d0*h*h+a*a))     ! From original Wait paper
   
  integrand=(exp_term/(u1+u2))*cos(l*a)   ! From original Wait paper
  QmjP=QmjP+integrand*dlambda(i)
  
  integrand2=k1*k1*(exp_term/(k2*k2*u1+k1*k1*u2))*cos(l*a)   ! From Olsen paper
  NmjM=NmjM+integrand2*dlambda(i)

! Wait, equation 17 integration  
  exp_term=exp(-2d0*u1*h)
  integrand3=2d0*( (l*l-u1*u2)/(k2*k2*u1+k1*k1*u2) )*exp_term*cos(l*a)     ! note sqrt removed as I think it is an error
  Integral_term=Integral_term+integrand3*dlambda(i)
  
  if (plot_integrand) then
    write(20,'(5ES12.4)')real(l),real(integrand),imag(integrand),real(integrand2),imag(integrand2)
  end if
  
end do

if (plot_integrand) then
  write(*,*)'Wait integrals, f=',real(f),' ALPHA=',cmplx(ALPHA)
  write(*,*)'Carson Approximation: ALPHA~',cmplx(log(2d0*h/a))
  write(*,*)' QmjP =',cmplx(QmjP),' NmjM=',cmplx(NmjM)
  write(*,*)'Carson~',cmplx(Jc/2d0),' NmjM~',cmplx(0d0)
  
  write(*,*)'sqrt(k1*k1-beta*beta)=',sqrt(k1*k1-beta*beta)
  
  close(unit=20)
  
  if(show_plots) then
    command='gnuplot i2plot.plt'
    CALL execute_command_line(command,EXITSTAT=exit_stat,CMDSTAT=cmd_stat)
  end if
  
!  STOP 1
  
end if

END SUBROUTINE Wait_integral2
