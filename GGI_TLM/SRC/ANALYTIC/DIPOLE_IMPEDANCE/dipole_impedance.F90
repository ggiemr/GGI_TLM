PROGRAM dipole

! THIS NEEDS CHECKING - RESULTS DON'T LOOK GOOD AT HIGH FREQUENCY

IMPLICIT NONE

real l,a,ftest

real f,fmin,fmax,fstep,lambda

real s,smin,smax,sstep

complex Z

integer nf,ns,i

! START

l=0.028
a=0.001

ftest=5.357E9
!ftest=4.0E9

nf=2500
ns=200

fmin=0.1e9
fmax=20e9
fstep=(fmax-fmin)/real(nf)

smin=0.001
smax=0.100
sstep=(smax-smin)/real(ns)

! self impedance as a function of frequency 

open(unit=10,file='Zs.dat')

! loop over frequency
do i=1,nf

  f=fmin+(i-1)*fstep
  lambda=2.998e8/f
  
!  a=lambda*0.001/2.0

  CALL dipole_self_impedance(l,a,f,Z)
  
!  write(10,*)l/lambda,real(Z),imag(Z)
  write(10,*)f,real(Z),imag(Z)
  
end do

close(unit=10)


! mutual impedance as a function of separation

open(unit=10,file='Zm.dat')

! loop over separation
do i=1,ns

  s=smin+(i-1)*sstep
  f=ftest
  
  CALL dipole_mutual_impedance(l,s,f,Z)
  
  write(10,*)s,real(Z),imag(Z)

end do

close(unit=10)

END PROGRAM dipole
!
! ____________________________________________________________
!
!
SUBROUTINE dipole_self_impedance(l,a,f,Z)

! Uses formula from Brown and King paper

IMPLICIT NONE

real    :: l,a,f
complex :: Z

real    :: k0
real    :: G
real*8  :: x1,ci1,si1
real*8  :: x2,ci2,si2
real*8  :: sc
real    :: cot
real*8  :: R,X

real,parameter :: gamma=0.5772156649    ! Euler constant

! START

k0=6.283185*f/2.998E8

! arguments of sine and cosine integrals
G=k0*l/2.0

x1=2.0*G
x2=4.0*G

CALL cisia ( x1, ci1, si1 )
CALL cisia ( x2, ci2, si2 )

sc=(sin(k0*l/2.0))**2
cot=cos(G)/sin(G)

R=30.0*(  (1.0-cot*cot)*( gamma+log(4.0*G)-ci2 )   &
         + 4.0*cot*cot* (gamma+log(2.0*G)-Ci1)     &
         + 2.0*cot*     (si2-2.0*si1)  )

X=30.0*(  2.0*cot*( gamma-log(l/(2.0*k0*a*a))+Ci2-2.0*ci1 ) &
         -(cot*cot-1.0)*(si2-2.0*si1)+2.0*si1/sc  )
         
Z=cmplx(R,X)

END SUBROUTINE dipole_self_impedance
!
! ____________________________________________________________
!
!
SUBROUTINE dipole_self_impedance_wikipedia(l,a,f,Z)

! THIS NEEDS CHECKING - RESULTS DON'T agree with the Brown and King data.

IMPLICIT NONE

real    :: l,a,f
complex :: Z

real    :: k0
real*8  :: x1,ci1,si1
real*8  :: x2,ci2,si2
real*8  :: x3,ci3,si3
real*8  :: sc
real*8  :: R,X

real,parameter :: gamma=0.5772156649    ! Euler constant

! START

k0=6.283185*f/2.998E8

! arguments of sine and cosine integrals

x1=k0*l
x2=2.0*k0*l
x3=2.0*k0*a*a/l

CALL cisia ( x1, ci1, si1 )
CALL cisia ( x2, ci2, si2 )
CALL cisia ( x3, ci3, si3 )

sc=(sin(k0*l/2.0))**2

R=(1.0/(6.283185*8.85e-12*2.998E8*sc))*( gamma+log(x1)-ci1+0.5*sin(x1)*(si2-2.0*si1)  &
                                        +0.5*cos(x1)*(gamma+log(x1/2.0)+ci2-2.0*ci1) )

X=(1.0/(2.0*6.283185*8.85e-12*2.998E8*sc))*( 2.0*si1+cos(x1)*(2.0*si1-si2)   & 
                                            -sin(x1)*(2.0*ci1-ci2-2.0*ci3) )  ! check sign here...

Z=cmplx(R,X)

END SUBROUTINE dipole_self_impedance_wikipedia
!
! ____________________________________________________________
!
!
SUBROUTINE dipole_mutual_impedance(l,d,f,Z)

! Brown and king

IMPLICIT NONE

real    :: l,d,f
complex :: Z
real    :: a,G

real    :: k0
real*8  :: x1,ci1,si1
real*8  :: x2,ci2,si2
real*8  :: x3,ci3,si3
real*8  :: x4,ci4,si4
real*8  :: x5,ci5,si5

real*8  :: R,X

! START

k0=6.283185*f/2.998E8

! arguments of sine and cosine integrals

a=l/2.0
G=k0*a

x1=k0*d
x2=k0*(sqrt(d*d+a*a)-a)
x3=k0*(sqrt(d*d+a*a)+a)
x4=k0*(sqrt(d*d+4.0*a*a)-2.0*a)
x5=k0*(sqrt(d*d+4.0*a*a)+2.0*a)

CALL cisia ( x1, ci1, si1 )
CALL cisia ( x2, ci2, si2 )
CALL cisia ( x3, ci3, si3 )
CALL cisia ( x4, ci4, si4 )
CALL cisia ( x5, ci5, si5 )

R=(30.0/(sin(G)**2))*(  2.0*( 2.0+cos(2.0*G))*ci1-4.0*cos(G)*cos(G)*( Ci2   &
                       +Ci3)+cos(2.0*G)*( Ci4      &
                       +Ci5 ) +sin(2.0*G)*( Si5   &
                       -Si4-2.0*si3               &
                       +2.0*si2 ) )

X=(30.0/(sin(G)**2))*( -2.0*( 2.0+cos(2.0*G))*si1+4.0*cos(G)*cos(G)*( Si2   &
                       +Si3 )-cos(2.0*G)*( Si4      &
                       +Si5 ) +sin(2.0*G)*( Ci5   &
                       -Ci4-2.0*Ci3               &
                       +2.0*Ci2 ) )


Z=cmplx(R,X)

END SUBROUTINE dipole_mutual_impedance
!
! ____________________________________________________________
!
!
SUBROUTINE dipole_mutual_impedance_sendy(l,s,f,Z)

! THIS NEEDS CHECKING - RESULTS DON'T LOOK GOOD AT HIGH FREQUENCY

IMPLICIT NONE

real    :: l,s,f
complex :: Z


real    :: k0
real*8  :: x1,ci1,si1
real*8  :: x2,ci2,si2
real*8  :: x3,ci3,si3

real*8  :: R,X

! START

k0=6.283185*f/2.998E8

! arguments of sine and cosine integrals

x1=k0*s
x2=k0*(sqrt(s*s+l*l)+l)
x3=k0*(sqrt(s*s+l*l)-l)

CALL cisia ( x1, ci1, si1 )
CALL cisia ( x2, ci2, si2 )
CALL cisia ( x3, ci3, si3 )

R=(1.0/(2.0*6.283185*8.85e-12*2.998E8))*(2.0*ci1-ci2-ci3)

X=-(1.0/(2.0*6.283185*8.85e-12*2.998E8))*(2.0*si1-si2-si3)

Z=cmplx(R,X)

END SUBROUTINE dipole_mutual_impedance_sendy

subroutine cisia ( x, ci, si )

!*****************************************************************************80
!
!! CISIA computes cosine Ci(x) and sine integrals Si(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Ci(x) and Si(x).
!
!    Output, real ( kind = 8 ) CI, SI, the values of Ci(x) and Si(x).
!
  implicit none

  real ( kind = 8 ) bj(101)
  real ( kind = 8 ) ci
  real ( kind = 8 ) el
  real ( kind = 8 ) eps
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) p2
  real ( kind = 8 ) si
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) xa
  real ( kind = 8 ) xa0
  real ( kind = 8 ) xa1
  real ( kind = 8 ) xcs
  real ( kind = 8 ) xf
  real ( kind = 8 ) xg
  real ( kind = 8 ) xg1
  real ( kind = 8 ) xg2
  real ( kind = 8 ) xr
  real ( kind = 8 ) xs
  real ( kind = 8 ) xss

  p2 = 1.570796326794897D+00
  el = 0.5772156649015329D+00
  eps = 1.0D-15
  x2 = x * x

  if ( x == 0.0D+00 ) then

    ci = -1.0D+300
    si = 0.0D+00

  else if ( x <= 16.0D+00 ) then

    xr = -0.25D+00 * x2
    ci = el + log ( x ) + xr
    do k = 2, 40
      xr = -0.5D+00 * xr * ( k - 1 ) / ( k * k * ( 2 * k - 1 ) ) * x2
      ci = ci + xr
      if ( abs ( xr ) < abs ( ci ) * eps ) then
        exit
      end if
    end do

    xr = x
    si = x
    do k = 1, 40
      xr = -0.5D+00 * xr * ( 2 * k - 1 ) / k / ( 4 * k * k + 4 * k + 1 ) * x2
      si = si + xr
      if ( abs ( xr ) < abs ( si ) * eps ) then
        return
      end if
    end do

  else if ( x <= 32.0D+00 ) then

    m = int ( 47.2D+00 + 0.82D+00 * x )
    xa1 = 0.0D+00
    xa0 = 1.0D-100
    do k = m, 1, -1
      xa = 4.0D+00 * k * xa0 / x - xa1
      bj(k) = xa
      xa1 = xa0
      xa0 = xa
    end do
    xs = bj(1)
    do k = 3, m, 2
      xs = xs + 2.0D+00 * bj(k)
    end do
    bj(1) = bj(1) / xs
    do k = 2, m
      bj(k) = bj(k) / xs
    end do
    xr = 1.0D+00
    xg1 = bj(1)
    do k = 2, m
      xr = 0.25D+00 * xr * ( 2.0D+00 * k - 3.0D+00 ) **2 &
        / ( ( k - 1.0D+00 ) * ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) * x
      xg1 = xg1 + bj(k) * xr
    end do

    xr = 1.0D+00
    xg2 = bj(1)
    do k = 2, m
      xr = 0.25D+00 * xr * ( 2.0D+00 * k - 5.0D+00 )**2 &
        / ( ( k-1.0D+00 ) * ( 2.0D+00 * k - 3.0D+00 ) ** 2 ) * x
      xg2 = xg2 + bj(k) * xr
    end do

    xcs = cos ( x / 2.0D+00 )
    xss = sin ( x / 2.0D+00 )
    ci = el + log ( x ) - x * xss * xg1 + 2.0 * xcs * xg2 - 2.0 * xcs * xcs
    si = x * xcs * xg1 + 2.0 * xss * xg2 - sin ( x )

  else

    xr = 1.0D+00
    xf = 1.0D+00
    do k = 1, 9
      xr = -2.0D+00 * xr * k * ( 2 * k - 1 ) / x2
      xf = xf + xr
    end do
    xr = 1.0D+00 / x
    xg = xr
    do k = 1, 8
      xr = -2.0D+00 * xr * ( 2 * k + 1 ) * k / x2
      xg = xg + xr
    end do
    ci = xf * sin ( x ) / x - xg * cos ( x ) / x
    si = p2 - xf * cos ( x ) / x - xg * sin ( x ) / x

  end if

  return
end

subroutine cisib ( x, ci, si )

!*****************************************************************************80
!
!! CISIB computes cosine and sine integrals.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    20 March 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Ci(x) and Si(x).
!
!    Output, real ( kind = 8 ) CI, SI, the values of Ci(x) and Si(x).
!
  implicit none

  real ( kind = 8 ) ci
  real ( kind = 8 ) fx
  real ( kind = 8 ) gx
  real ( kind = 8 ) si
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  x2 = x * x

  if ( x == 0.0D+00 ) then

    ci = -1.0D+300
    si = 0.0D+00

  else if ( x <= 1.0D+00 ) then

    ci = (((( -3.0D-08        * x2 &
             + 3.10D-06     ) * x2 &
             - 2.3148D-04   ) * x2 &
             + 1.041667D-02 ) * x2 &
             - 0.25D+00     ) * x2 + 0.577215665D+00 + log ( x )

     si = (((( 3.1D-07        * x2 &
             - 2.834D-05    ) * x2 &
             + 1.66667D-03  ) * x2 &
             - 5.555556D-02 ) * x2 + 1.0D+00 ) * x

  else

    fx = (((( x2              &
      + 38.027264D+00  ) * x2 &
      + 265.187033D+00 ) * x2 &
      + 335.67732D+00  ) * x2 &
      + 38.102495D+00  ) /    &
      (((( x2                 &
      + 40.021433D+00  ) * x2 &
      + 322.624911D+00 ) * x2 &
      + 570.23628D+00  ) * x2 &
      + 157.105423D+00 )

    gx = (((( x2               &
      + 42.242855D+00  ) * x2  &
      + 302.757865D+00 ) * x2  &
      + 352.018498D+00 ) * x2  &
      + 21.821899D+00 ) /      &
      (((( x2                  &
      + 48.196927D+00   ) * x2 &
      + 482.485984D+00  ) * x2 &
      + 1114.978885D+00 ) * x2 &
      + 449.690326D+00  ) / x

    ci = fx * sin ( x ) / x - gx * cos ( x ) / x

    si = 1.570796327D+00 - fx * cos ( x ) / x - gx * sin ( x ) / x

  end if

  return
end
