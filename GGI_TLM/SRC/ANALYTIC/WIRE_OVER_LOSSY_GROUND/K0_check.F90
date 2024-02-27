! Check the integral form of the K0 function

open(unit=20,file='expansion_test.dat')

do ix=1000,3000

! point
  xp=dble(ix)/100d0
  yp=0.1
  rp=sqrt((xp-h)*(xp-h)+yp*yp)

  beta=(3.758723289E-02,-5.500228726E-04)

  kr=sqrt(k1*k1-beta*beta)

! get the discretisation for l (lambda in Wait paper)
  
  nlambda_local=nlambda*2

  fimin=1d-4
  lmax=-20.0*log(fimin)/(2d0*h)
  dl=2d0*lmax/dble(nlambda_local)

  write(*,*)'lmax=',lmax,'dl=',dl

  kr=sqrt(k1*k1-beta*beta)
  if(imag(kr).GT.0d0) kr=kr*(-1d0)   ! sign check for lossy propagation away from the wire

  arg=rp*j*kr ! argument of modified Bessel function
  CALL calc_modified_Bessel_function((0d0,0d0),arg,K0ra)

! Integral in Wait paper, equation 4

  K0r=(0d0,0d0)

  do i=-nlambda_local,nlambda_local

    l=dble(i)*dl

    u1=sqrt(l*l+beta*beta-k1*k1)
    if(dble(u1).LT.0d0) u1=u1*(-1d0)   ! sign check. Olsen paper section III
  
    if (xp.lt.h) then
      K0r=K0r+( 0.5d0*exp(-j*l*yp)*exp( u1*(xp-h))/u1 )*dl
    else
      K0r=K0r+( 0.5d0*exp(-j*l*yp)*exp(-u1*(xp-h))/u1 )*dl
    end if
    
  end do

!write(*,*)'K0r=',K0r

  write(20,'(5ES12.4)')xp,real(k0ra),imag(K0ra),real(k0r),imag(K0r)

end do ! next x position

close(unit=20)

STOP 1

K0r=(0d0,0d0)

do i=1,nlambda

  l=dble(i)*dl-dl/2d0

  u1=sqrt(l*l+beta*beta-k1*k1)
  if(dble(u1).LT.0d0) u1=u1*(-1d0)   ! sign check. Olsen paper section III
  
  K0r=K0r+( cos(rp*l)/u1 )*dl
  
!  write(20,'(3ES12.4)')real(l),exp(-j*l*yp)*exp(u1*(xp-h))/u1*dl
!  write(20,'(3ES12.4)')real(l),exp(-j*l*yp)
  
end do

write(*,*)'CJS cosine form:'
write(*,*)'K0r=',K0r

K0r=(0d0,0d0)

do i=-nlambda,nlambda

  l=dble(i)*dl

  u1=sqrt(l*l+beta*beta-k1*k1)
  if(dble(u1).LT.0d0) u1=u1*(-1d0)   ! sign check. Olsen paper section III
  
  K0r=K0r+0.5d0*( exp(j*rp*l)/u1 )*dl
  
!  write(20,'(3ES12.4)')real(l),exp(-j*l*yp)*exp(u1*(xp-h))/u1*dl
!  write(20,'(3ES12.4)')real(l),exp(-j*l*yp)
 
end do

write(*,*)'CJS exponential form:'
write(*,*)'K0r=',K0r

close(unit=20)



STOP 1
