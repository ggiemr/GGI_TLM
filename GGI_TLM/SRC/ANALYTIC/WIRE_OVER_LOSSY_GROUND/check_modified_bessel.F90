
write(*,*)'Check implementation of modified Bessel functions'

kr=sqrt(k1*k1-beta*beta)
if(imag(kr).GT.0d0) kr=kr*(-1d0)   ! sign check for lossy propagation away from the wire

xp=mcx 
yp=mcy 

write(*,*)'Point:',real(xp),real(yp)


! Analytic calculation of modified bessel function and partial derivatives wrt x and y

r=sqrt( (xp-h)**2+yp**2)

write(*,*)'r=',real(r)

write(*,*)'kr=',cmplx(kr)

arg=(0d0,1d0)*r*kr ! argument of modified Bessel function
write(*,*)'arg=',cmplx(arg)
CALL calc_modified_Bessel_function_and_derivative((0d0,0d0),arg,K0r,dK0r)

partial_K0r_x=dK0r*(0d0,1d0)*kr*(xp-h)/r
partial_K0r_y=dK0r*(0d0,1d0)*kr*yp/r

! Calculate modified Bessel function and partial derivatives from the Integral form
K0r_i          =(0d0,0d0)
partial_K0r_x_i=(0d0,0d0)
partial_K0r_y_i=(0d0,0d0)

dl=2d0*lmax/dble(nlambda)
write(*,*)'lmax=',lmax,'dl=',dl

  K0r=(0d0,0d0)
  do i=-nlambda,nlambda
    l=dble(i)*dl
    u1=sqrt(l*l+beta*beta-k1*k1)
    if(dble(u1).LT.0d0) u1=u1*(-1d0)   ! sign check. Olsen paper section III  
    if (xp.lt.h) then
      K0r=K0r+( 0.5d0*exp(-j*l*yp)*exp( u1*(xp-h))/u1 )*dl
    else
      K0r=K0r+( 0.5d0*exp(-j*l*yp)*exp(-u1*(xp-h))/u1 )*dl
    end if    
  end do

write(*,*)'integral1, K0r=',K0r

do i=-nlambda,nlambda

!  l=lambda_field(i)
!  dl=dl_integral_field(i)

  l=dble(i)*dl
  u1=sqrt(l*l+beta*beta-k1*k1)
  if(dble(u1).LT.0d0) u1=u1*(-1d0)   ! sign check. Olsen paper section III

  if (xp.lt.h) then
    Pe_px_my=0.5d0*( exp( u1*(xp-h)) )*(exp(-j*l*yp)/u1)*dl  
    K0r_i          =K0r_i           + Pe_px_my   
    partial_K0r_x_i=partial_K0r_x_i + Pe_px_my*u1
    partial_K0r_y_i=partial_K0r_y_i + Pe_px_my*(-j*l)
  else
    Pe_px_my=0.5d0*( exp(-u1*(xp-h)) )*(exp(-j*l*yp)/u1)*dl  
    K0r_i          =K0r_i           + Pe_px_my   
    partial_K0r_x_i=partial_K0r_x_i - Pe_px_my*u1
    partial_K0r_y_i=partial_K0r_y_i + Pe_px_my*(-j*l)
  end if
  
        
end do ! next l (lambda)

write(*,*)'mod_zbes : K0r=',K0r
write(*,*)'integral : K0r=',K0r_i

write(*,*)'mod_zbes : dK0r/dx=',partial_K0r_x
write(*,*)'integral : dK0r/dx=',partial_K0r_x_i

! numerical differential
delta_l=0.00001d0

xp=mcx+delta_l/2d0
yp=mcy 
r=sqrt( (xp-h)**2+yp**2)
arg=(0d0,1d0)*r*kr ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0rp)

xp=mcx-delta_l/2d0
yp=mcy 
r=sqrt( (xp-h)**2+yp**2)
arg=(0d0,1d0)*r*kr ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0rm)

write(*,*)'numerical: dK0r/dx=',(K0rp-K0rm)/delta_l


write(*,*)'mod_zbes : dK0r/dy=',partial_K0r_y
write(*,*)'integral : dK0r/dy=',partial_K0r_y_i

xp=mcx
yp=mcy+delta_l/2d0
r=sqrt( (xp-h)**2+yp**2)
arg=(0d0,1d0)*r*kr ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0rp)

xp=mcx
yp=mcy-delta_l/2d0
r=sqrt( (xp-h)**2+yp**2)
arg=(0d0,1d0)*r*kr ! argument of modified Bessel function
CALL calc_modified_Bessel_function((0d0,0d0),arg,K0rm)

write(*,*)'numerical: dK0r/dy=',(K0rp-K0rm)/delta_l

STOP 1
