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
!
!_____________________________________________________
!       
       
       subroutine rcs(theta,phi,beta,a,nmax,sigma3d)

USE constants

       implicit none
              
       real*8 theta,phi,beta,a
       real*8 sigma3d
       integer nmax
       
! bessel functions       
       complex*16 sp_besselJ(0:nmax)
       complex*16 sp_bessely(0:nmax)
       complex*16 sp_besselH1(0:nmax)
       complex*16 sp_besselH2(0:nmax)
       complex*16 sp_dbesselJ(0:nmax)
       complex*16 sp_dbessely(0:nmax)
       complex*16 sp_dbesselH1(0:nmax)
       complex*16 sp_dbesselH2(0:nmax)
       
! Legendre polynomials
       real*8 Pn0,Pn1,dPn1
       
       complex*16 betar
       
       complex*16 an,bn,cn
       complex*16 Etheta,Ephi
       
       real*8 ct,st
       
       integer n

! START
       
! set parameters       
       betar=dcmplx(beta*a)

! reset fields
       Etheta=(0D0,0D0)
       Ephi  =(0D0,0D0)
       
! Evaluate bessel functions    

       CALL spherical_bessel2( betar, nmax, sp_besselJ, sp_dbesselJ, sp_bessely, sp_dbessely )
       
! normal definition of spherical Bessel/ Hankel functions
       sp_besselH1(0:nmax)=sp_besselJ(0:nmax)+j*sp_bessely(0:nmax)
       sp_besselH2(0:nmax)=sp_besselJ(0:nmax)-j*sp_bessely(0:nmax)
       
       sp_dbesselH1(0:nmax)=sp_dbesselJ(0:nmax)+j*sp_dbessely(0:nmax)
       sp_dbesselH2(0:nmax)=sp_dbesselJ(0:nmax)-j*sp_dbessely(0:nmax)
       
       do n=1,nmax
       
         an=(j**(-n))*(2d0*n+1)/(n*(n+1))
	 
	 bn=-an*sp_dbesselJ(n)/sp_dbesselH2(n)
     
	 cn=-an*sp_besselJ(n)/sp_besselH2(n)
     
! cope with the limiting case of theta=0,180 degrees...     
	 ct=cos(theta)	 
	 st=sin(theta)
	 if(ct.eq. 1D0)  ct= 1D0-0.00001D0
	 if(ct.eq.-1D0) ct=-1D0+0.00001D0
         st=sqrt(1d0-ct*ct)
              
         Etheta=Etheta+j*(exp(-j*betar)/beta)*cos(phi)*(j**n)		&
                        *( bn*st*dPn1(n,ct) -cn*Pn1(n,ct)/st  )
     
         Ephi=Ephi+j*( exp(-j*betar)/beta )*sin(phi)*(j**n)		&
                    *(  bn*Pn1(n,ct)/st -cn*st*dPn1(n,ct)  )
	 
       end do
       
       sigma3d=   4D0*pi*( abs(Etheta)*abs(Etheta) + abs(Ephi)*abs(Ephi)  )
              
       return       
       end
