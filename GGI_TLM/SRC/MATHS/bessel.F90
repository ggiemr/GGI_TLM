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
! Evaluation of sphereical bessel functions, alternative spherical bessel functions and
! functions for derivatives of bessel funcrtions.

! The bessel function evaluation is done using Fortran90 intrinsic functions

! The spherical bessel function evaluation is done using a simple recurrence relation. This may
! not be robust and should be used with care. Ideally this should be updated with a more robust
! subroutine. 

SUBROUTINE spherical_bessel( z, nmax, Jn, dJn, Yn, dYn )

  integer 	:: nmax
  complex*16 	:: z
  
  complex*16	:: Jn(0:nmax),dJn(0:nmax),Yn(0:nmax),dYn(0:nmax)

! local variables

  integer 	:: n_starting_iterations
  complex*16	:: Jn_temp_pp,Jn_temp_p,Jn_temp,J0,Jscale
    
  complex*16,parameter	:: j=(0d0,1d0)
       
  integer n,i,constantJ,constantY
    
! START
       
!  if (abs(z).lt.1d-40) then     
!      
!    do n=0,2
!  
!      constantJ=1
!      constantY=-1
!      do i=1,n
!        constantJ=constantJ*(2*i+1)
!        constantY=constantY*(2*i-1)
!      end do
!  
!      Jn(n)=(z**n)/constantJ
!    
!      Yn(n)=constantY/(z**(n+1))
!      
!      dJn(n)=n*(z**(n-1))/constantJ   
!      dYn(n)=-(n+1)*constantY/(z**(n+2))
!    
!    end do
!    
!  else

! Evaluate functions from recurrence relations  
! direct evaluation of order 0 and 1 functions
    J0 = sin(z)/z
  
! backward recurrence for Jn functions
    n_starting_iterations=int(abs(z))+nmax+2   
    Jn_temp_pp= 1d-100
    Jn_temp_p = 1d-99

! do a few iterations of the reverse recurrence to get some starting values  
    do n=nmax+n_starting_iterations,nmax-1,-1
      Jn_temp=((2d0*n+3d0)/z)*Jn_temp_p-Jn_temp_pp
      Jn_temp_pp=Jn_temp_p
      Jn_temp_p=Jn_temp      
    end do
  
    Jn(nmax)=(Jn_temp_pp/Jn_temp)*1d-100
    Jn(nmax-1)=(Jn_temp_p/Jn_temp)*1d-100
    
    do n=nmax-2,0,-1
      Jn(n)=((2d0*n+3d0)/z)*Jn(n+1)-Jn(n+2)
    end do

! scale all terms in the series to give the known J0 and hence all the other values correctly
    Jscale=J0/Jn(0)
    Jn(0:nmax)=Jn(0:nmax)*Jscale
  
    Yn(0) = -cos(z)/z	 
    Yn(1) = (Yn(0)-sin(z))/z

! forward recurrence for Yn 
    do n=2,nmax
  
      Yn(n) = ((2d0*n-1)/z)*Yn(n-1)-Yn(n-2)
  
    end do ! next n
  
    dJn(0)=( cos(z)-( sin(z)/z ) )/z
    dJn(1)=Jn(0)-(2d0/z)*Jn(1)
    dYn(0)=( sin(z)+( cos(z)/z ) )/z
    dYn(1)=Yn(0)-(2d0/z)*Yn(1)

! forward recurrence for derivative functions  
    do n=2,nmax

      dJn(n)= Jn(n-1)-((n+1)/z)*Jn(n)
      dYn(n)= Yn(n-1)-((n+1)/z)*Yn(n)
  
    end do ! next n
  
!  end if ! not small argument

END SUBROUTINE spherical_bessel
!
! _________________________________________________________________
!

! Following is the evaluation of alternative spherical bessel functions
! These are used in the calculation of RCS of spheres. 
!
! The alternative spherical bessel functions are given by Fn(z)=zfn(z)

SUBROUTINE spherical_bessel2( z, nmax, Jn, dJn, Yn, dYn )

  integer 	:: nmax
  complex*16 	:: z
  
  complex*16	:: Jn(0:nmax),dJn(0:nmax),Yn(0:nmax),dYn(0:nmax)

! local variables

  integer 	:: n_starting_iterations
  complex*16	:: Jn_temp_pp,Jn_temp_p,Jn_temp,J0,Jscale
    
  complex*16,parameter	:: j=(0d0,1d0)
    
! START

  CALL spherical_bessel( z, nmax, Jn, dJn, Yn, dYn )

! return a function F(z)=zf(z)  dF(z)/dz= zdf(z)/dz+f(z)  
  do n=0,nmax

    dJn(n)=z*dJn(n)+Jn(n)
    dYn(n)=z*dYn(n)+Yn(n)
  
    Jn(n) =z*Jn(n)
    Yn(n) =z*Yn(n)
  
  end do ! next n

END SUBROUTINE spherical_bessel2
!
! _________________________________________________________________
!

FUNCTION derivative_besselJ(n,x) result(res)

integer n
real*8 x
real*8 res

! START

  res=-besjn(n+1,x)+(n/x)*besjn(n,x)

END FUNCTION derivative_besselJ
!
! _________________________________________________________________
!

FUNCTION derivative_besselY(n,x) result(res)

integer n
real*8 x
real*8 res

! START

  res=-besyn(n+1,x)+(n/x)*besyn(n,x)

END FUNCTION derivative_besselY

! 

