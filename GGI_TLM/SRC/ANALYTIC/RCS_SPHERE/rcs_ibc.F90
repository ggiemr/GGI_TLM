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
       
       subroutine rcs_ibc(theta,phi,w,a0,nmax,z11,sigma3d)

 USE Constants

       implicit none
              
       real*8 theta,phi,w,a0
       complex*16 z11
       integer nmax
       complex*16 Etheta,Ephi,Htheta,Hphi
       real*8 sigma3d
     
! bessel functions       
       
       complex*16 sp_besselJ_beta0_a0(0:nmax)
       complex*16 sp_bessely_beta0_a0(0:nmax)
       complex*16 sp_besselH1_beta0_a0(0:nmax)
       complex*16 sp_besselH2_beta0_a0(0:nmax)
       complex*16 sp_dbesselJ_beta0_a0(0:nmax)
       complex*16 sp_dbessely_beta0_a0(0:nmax)
       complex*16 sp_dbesselH1_beta0_a0(0:nmax)
       complex*16 sp_dbesselH2_beta0_a0(0:nmax)
       
! Legendre polynomials
       real*8 Pn0,Pn1,dPn1
       
! constants       
       real*8 lambda
       
       complex*16 beta0
       complex*16 beta0_a0
       
       complex*16 z2
       
! matrix elements
       complex*16 k11,k1,k22,k2     
              
       complex*16 an,bn,cn
       
       complex*16 Ethetai,Ephii,Eri
       complex*16 Ethetaf,Ephif

       real*8 ct,st
              
       integer n
       
       real*8 oterm,iterm,osum,isum,norm
       
! set parameters       

       beta0=cmplx(w/c0)
              
       beta0_a0=beta0*a0
       
! evaluate bessel functions

       CALL spherical_bessel2( beta0_a0, nmax, sp_besselJ_beta0_a0,	&
                                              sp_dbesselJ_beta0_a0,	&
					      sp_bessely_beta0_a0,	&
					      sp_dbessely_beta0_a0 )
					   
       sp_besselH1_beta0_a0(0:nmax)=sp_besselJ_beta0_a0(0:nmax)+j*sp_bessely_beta0_a0(0:nmax)
       sp_besselH2_beta0_a0(0:nmax)=sp_besselJ_beta0_a0(0:nmax)-j*sp_bessely_beta0_a0(0:nmax)
       sp_dbesselH1_beta0_a0(0:nmax)=sp_dbesselJ_beta0_a0(0:nmax)+j*sp_dbessely_beta0_a0(0:nmax)
       sp_dbesselH2_beta0_a0(0:nmax)=sp_dbesselJ_beta0_a0(0:nmax)-j*sp_dbessely_beta0_a0(0:nmax)
              
! reset fields
       Etheta=(0.0,0.0)
       Ephi=(0.0,0.0)
       Htheta=(0.0,0.0)
       Hphi=(0.0,0.0)

       Ethetai=(0.0,0.0)
       Ephii=(0.0,0.0)
       
       Ethetaf=(0.0,0.0)
       Ephif=(0.0,0.0)
       
!       small=1d-6  ! set in constants
       osum=0.0
       isum=0.0
       
       do n=1,nmax

! incident field expansion coefficient       
         an=(j**(-n))*(2.0*n+1)/(n*(n+1))
	  
	 k11=(-beta0/(j*w*w*mu0*eps0))*sp_dbesselH2_beta0_a0(n)	&
            -z11*(-1.0/(w*mu0))*sp_besselH2_beta0_a0(n)
     
	 k1=-((-beta0/(j*w*w*mu0*eps0))*sp_dbesselJ_beta0_a0(n)	&
            -z11*(-1.0/(w*mu0))*sp_besselJ_beta0_a0(n) )
	 
	 k22=(1.0/(w*eps0*z0))*sp_besselH2_beta0_a0(n)	&
            -z11*(beta0/(j*w*w*mu0*eps0*z0))*sp_dbesselH2_beta0_a0(n)
     
	 k2=-( (1.0/(w*eps0*z0))*sp_besselJ_beta0_a0(n)	&
            -z11*(beta0/(j*w*w*mu0*eps0*z0))*sp_dbesselJ_beta0_a0(n) ) 
	 
! cope with the limiting case of theta=0,180 degrees...     
	 ct=cos(theta)	 
	 if(ct.eq.1.0) ct=1.0-0.0000001
	 if(ct.eq.-1.0) ct=-1.0+0.0000001
         st=sqrt(1.0-ct*ct)
     
! given the expansion coefficients, calculate the scattered E field in the far field    
           Ethetaf=Ethetaf+j*(exp(-j*beta0_a0)/beta0)*cos(phi)*(j**n)	&
                *(  bn*st*dPn1(n,ct)	&
                   -cn*Pn1(n,ct)/st    )
     
           Ephif=Ephif+j*(exp(-j*beta0_a0)/beta0)*sin(phi)*(j**n)	&
                *(  bn*Pn1(n,ct)/st	&
                   -cn*st*dPn1(n,ct)   )
     

! Convergence testing
         oterm=sqrt(abs(bn*bn)+abs(cn*cn))
	 osum=osum+oterm
	 
	 if ((oterm)/(osum).lt.small) goto 1234
     
       end do
       
1234   continue       
        
       print*,'Terminated series, n=',n
       
       sigma3d=   4.0*pi*( abs(Ethetaf) *abs(Ethetaf) + abs(Ephif)   *abs(Ephif)    )

       print*,abs(Ethetai)+abs(Ephii)

! TOTAL FIELD         
       Etheta=Etheta
       Ephi=Ephi
! SCATTERED FIELD       
!       Etheta=Etheta-Ethetai
!       Ephi=Ephi-Ephii
              
       return       
       end
