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
       
       subroutine rcs_dielectric_ibc(theta,phi,w,a2,eps2,mu2,	&
                  z11,z12,z21,z22,nmax,Etheta,Ephi,Htheta,Hphi,sigma3d)

USE constants

       implicit none
              
       real*8 theta,phi,w,a2
       complex*16 eps2,mu2
       complex*16 z11,z12,z21,z22
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
       
       complex*16 sp_besselJ_beta0_a2(0:nmax)
       complex*16 sp_bessely_beta0_a2(0:nmax)
       complex*16 sp_besselH1_beta0_a2(0:nmax)
       complex*16 sp_besselH2_beta0_a2(0:nmax)
       complex*16 sp_dbesselJ_beta0_a2(0:nmax)
       complex*16 sp_dbessely_beta0_a2(0:nmax)
       complex*16 sp_dbesselH1_beta0_a2(0:nmax)
       complex*16 sp_dbesselH2_beta0_a2(0:nmax)
       
       complex*16 sp_besselJ_beta2_a0(0:nmax)
       complex*16 sp_bessely_beta2_a0(0:nmax)
       complex*16 sp_besselH1_beta2_a0(0:nmax)
       complex*16 sp_besselH2_beta2_a0(0:nmax)
       complex*16 sp_dbesselJ_beta2_a0(0:nmax)
       complex*16 sp_dbessely_beta2_a0(0:nmax)
       complex*16 sp_dbesselH1_beta2_a0(0:nmax)
       complex*16 sp_dbesselH2_beta2_a0(0:nmax)
       
       complex*16 sp_besselJ_beta2_a2(0:nmax)
       complex*16 sp_bessely_beta2_a2(0:nmax)
       complex*16 sp_besselH1_beta2_a2(0:nmax)
       complex*16 sp_besselH2_beta2_a2(0:nmax)
       complex*16 sp_dbesselJ_beta2_a2(0:nmax)
       complex*16 sp_dbessely_beta2_a2(0:nmax)
       complex*16 sp_dbesselH1_beta2_a2(0:nmax)
       complex*16 sp_dbesselH2_beta2_a2(0:nmax)
       
! Legendre polynomials
       real*8 Pn0,Pn1,dPn1
       
! constants       
       real*8 lambda
       
       complex*16 beta0,beta2
       real*8 a0
       complex*16 beta0_a2,beta2_a2,beta2_a0,beta0_a0
       
       complex*16 z2
       
! matrix elements
       complex*16 k3b,k3d,k3e      
       complex*16 k4c,k4f,k4g      
       complex*16 k5b,k5d,k5e      
       complex*16 k6c,k6f,k6g
       complex*16 r3,r4,r5,r6
       
       complex*16 u11,u12,u21,u22,det,rdet   
       
       complex*16 an,bn,cn,dn,fn
       
       complex*16 Ethetai,Ephii,Eri
       complex*16 Ethetaf,Ephif

       real*8 ct,st
              
       integer n
       
       real*8 oterm,iterm,osum,isum,norm,cvg

       cvg=1d-6      ! convergence criterion

! small non zero radius for field calculation
       a0=a2/1000D0

! set parameters       

       beta0=cmplx(w/c0)
       beta2=w*sqrt(mu2*eps2)
       if (aimag(beta2).gt.0D0) beta2=-beta2
       z2=sqrt(mu2/eps2)
       if (real(beta2).lt.0D0) beta2=-beta2
              
       beta0_a2=beta0*a2
       beta0_a0=beta0*a0
       beta2_a2=beta2*a2
       beta2_a0=beta2*a0
       
! evaluate bessel functions

       CALL spherical_bessel2( beta0_a0, nmax, sp_besselJ_beta0_a0,	&
                                               sp_dbesselJ_beta0_a0,	&
					       sp_bessely_beta0_a0,	&
					       sp_dbessely_beta0_a0 )
					   
       sp_besselH1_beta0_a0(0:nmax)=sp_besselJ_beta0_a0(0:nmax)+j*sp_bessely_beta0_a0(0:nmax)
       sp_besselH2_beta0_a0(0:nmax)=sp_besselJ_beta0_a0(0:nmax)-j*sp_bessely_beta0_a0(0:nmax)
       sp_dbesselH1_beta0_a0(0:nmax)=sp_dbesselJ_beta0_a0(0:nmax)+j*sp_dbessely_beta0_a0(0:nmax)
       sp_dbesselH2_beta0_a0(0:nmax)=sp_dbesselJ_beta0_a0(0:nmax)-j*sp_dbessely_beta0_a0(0:nmax)

       CALL spherical_bessel2( beta0_a2, nmax, sp_besselJ_beta0_a2,	&
                                               sp_dbesselJ_beta0_a2,	&
					       sp_bessely_beta0_a2,	&
					       sp_dbessely_beta0_a2 )
					   
       sp_besselH1_beta0_a2(0:nmax)=sp_besselJ_beta0_a2(0:nmax)+j*sp_bessely_beta0_a2(0:nmax)
       sp_besselH2_beta0_a2(0:nmax)=sp_besselJ_beta0_a2(0:nmax)-j*sp_bessely_beta0_a2(0:nmax)
       sp_dbesselH1_beta0_a2(0:nmax)=sp_dbesselJ_beta0_a2(0:nmax)+j*sp_dbessely_beta0_a2(0:nmax)
       sp_dbesselH2_beta0_a2(0:nmax)=sp_dbesselJ_beta0_a2(0:nmax)-j*sp_dbessely_beta0_a2(0:nmax)

       CALL spherical_bessel2( beta2_a0, nmax, sp_besselJ_beta2_a0,	&
                                               sp_dbesselJ_beta2_a0,	&
					       sp_bessely_beta2_a0,	&
					       sp_dbessely_beta2_a0 )
					   
       sp_besselH1_beta2_a0(0:nmax)=sp_besselJ_beta2_a0(0:nmax)+j*sp_bessely_beta2_a0(0:nmax)
       sp_besselH2_beta2_a0(0:nmax)=sp_besselJ_beta2_a0(0:nmax)-j*sp_bessely_beta2_a0(0:nmax)
       sp_dbesselH1_beta2_a0(0:nmax)=sp_dbesselJ_beta2_a0(0:nmax)+j*sp_dbessely_beta2_a0(0:nmax)
       sp_dbesselH2_beta2_a0(0:nmax)=sp_dbesselJ_beta2_a0(0:nmax)-j*sp_dbessely_beta2_a0(0:nmax)

       CALL spherical_bessel2( beta2_a2, nmax, sp_besselJ_beta2_a2,	&
                                               sp_dbesselJ_beta2_a2,	&
					       sp_bessely_beta2_a2,	&
					       sp_dbessely_beta2_a2 )
					   
       sp_besselH1_beta2_a2(0:nmax)=sp_besselJ_beta2_a2(0:nmax)+j*sp_bessely_beta2_a2(0:nmax)
       sp_besselH2_beta2_a2(0:nmax)=sp_besselJ_beta2_a2(0:nmax)-j*sp_bessely_beta2_a2(0:nmax)
       sp_dbesselH1_beta2_a2(0:nmax)=sp_dbesselJ_beta2_a2(0:nmax)+j*sp_dbessely_beta2_a2(0:nmax)
       sp_dbesselH2_beta2_a2(0:nmax)=sp_dbesselJ_beta2_a2(0:nmax)-j*sp_dbessely_beta2_a2(0:nmax)
       
       lambda=2D0*pi*c0/w
       
! reset fields
       Etheta=(0D0,0D0)
       Ephi=(0D0,0D0)
       Htheta=(0D0,0D0)
       Hphi=(0D0,0D0)

       Ethetai=(0D0,0D0)
       Ephii=(0D0,0D0)
       
       Ethetaf=(0D0,0D0)
       Ephif=(0D0,0D0)
       
       osum=0D0
       isum=0D0
               
       
       do n=1,nmax

! incident field expansion coefficient       
         an=(j**(-n))*(2.0*n+1)/(n*(n+1))
	 
! Equations from Ephi, Htheta set
	 k3b=(-beta0/(j*w*w*mu0*eps0))*sp_dbesselH2_beta0_a2(n)	&
            -z11*(-1.0/(w*mu0))*sp_besselH2_beta0_a2(n)	
      
	 k3d=-z12*(-1.0/(w*mu2))*sp_besselJ_beta2_a2(n)
	 
	 k4c=(1.0/(w*eps0*z0))*sp_besselH2_beta0_a2(n)	&
       -z11*(beta0/(j*w*w*mu0*eps0*z0))*sp_dbesselH2_beta0_a2(n)	 
     
	 k4f=-z12*(beta2/(j*w*w*mu2*eps2*z2))*sp_dbesselJ_beta2_a2(n)
	 
	 k5b=-z21*(-1.0/(w*mu0))*sp_besselH2_beta0_a2(n)
	 
	 k5d=(-beta2/(j*w*w*mu2*eps2))*sp_dbesselJ_beta2_a2(n)	&
       -z22*(-1.0/(w*mu2))*sp_besselJ_beta2_a2(n)
	 
	 k6c=-z21*( beta0/(j*w*w*mu0*eps0*z0))*sp_dbesselH2_beta0_a2(n)	
	  
	 k6f= (+1.0/(w*eps2*z2))*sp_besselJ_beta2_a2(n)	&
     	 -z22*(beta2/(j*w*w*mu2*eps2*z2))*sp_dbesselJ_beta2_a2(n)	 

! RHS elements for solution of field expansion coefficients

         r3=-( an*(-beta0/(j*w*w*mu0*eps0))*sp_dbesselJ_beta0_a2(n)	&
         -z11*an*(-1.0/(w*mu0))*sp_besselJ_beta0_a2(n) )
     
	 r4=-( an*(1.0/(w*eps0*z0))*sp_besselJ_beta0_a2(n)	&
         -z11*an*(beta0/(j*w*w*mu0*eps0*z0))*sp_dbesselJ_beta0_a2(n) )
	 
	 r5=-( an*(-z21)*(-1.0/(w*mu0))*sp_besselJ_beta0_a2(n)  )
	 
	 r6=-( an*(-z21)*( beta0/(j*w*w*mu0*eps0*z0))*	&
               	   sp_dbesselJ_beta0_a2(n) )

! solution of submatrix equation 1
	
	 u11=k3b
	 u12=k3d
	 u21=k5b
	 u22=k5d
	 
	 norm=(abs(U11)+abs(U12)+abs(U21)+abs(U22))
	 
	 U11=U11/norm
	 U12=U12/norm
	 R3= R3 /norm
	 
	 U21=U21/norm
	 U22=U22/norm
	 R5= R5 /norm

         det=u11*u22-u12*u21
	 	 
	 if (abs(det).ne.0.0) then
	   rdet=1.0/det
	 else 
	   rdet=0.0
	 end if

	 bn=( u22*r3-u12*r5)*rdet
	 dn=(-u21*r3+u11*r5)*rdet
	 
! solution of submatrix equation 2	 
	
	 u11=k4c
	 u12=k4f
	 u21=k6c
	 u22=k6f
	 
	 norm=(abs(U11)+abs(U12)+abs(U21)+abs(U22))
	 
	 U11=U11/norm
	 U12=U12/norm
	 R4= R4 /norm
	 
	 U21=U21/norm
	 U22=U22/norm
	 R6= R6 /norm

         det=u11*u22-u12*u21
	 	 
	 if (abs(det).ne.0.0) then
	   rdet=1.0/det
	 else 
	   rdet=0.0
	 end if
	 cn=( u22*r4-u12*r6)*rdet
	 fn=(-u21*r4+u11*r6)*rdet
	 
! cope with the limiting case of theta=0,180 degrees...     
	 ct=cos(theta)	 
	 if(ct.eq.1.0) ct=1.0-0.0000001
	 if(ct.eq.-1.0) ct=-1.0+0.0000001
         st=sqrt(1.0-ct*ct)

! near field expression: total field in dielectric sphere

         Etheta=Etheta	&
        -dn* (beta2/(j*w*w*mu2*eps2*a0))*cos(phi)*	&
        sp_dbesselJ_beta2_a0(n)*dPn1(n,ct)*st-	&
        fn* ((1.0,0.0)/(w*eps2*z2*a0))*cos(phi)*	&
        sp_besselJ_beta2_a0(n)*Pn1(n,ct)/st     
     
         Ephi=Ephi+	&
        dn* (-beta2/(j*w*w*mu2*eps2*a0))*sin(phi)*	&
        sp_dbesselJ_beta2_a0(n)*Pn1(n,ct)/st-	&
        fn* ((1.0,0.0)/(w*eps2*z2*a0))*sin(phi)	&
        *sp_besselJ_beta2_a0(n)*dPn1(n,ct)*st     


! near field expression: incident field


         Ethetai=Ethetai	&
        -an* (beta0/(j*w*w*mu0*eps0*a0))*cos(phi)*	&
         sp_dbesselJ_beta0_a0(n)*dPn1(n,ct)*st-	&
        an* ((1.0,0.0)/(w*eps0*z0*a0))*cos(phi)*	&
        sp_besselJ_beta0_a0(n)*Pn1(n,ct)/st     

         Ephii=Ephii+	&
        an* (-beta0/(j*w*w*mu0*eps0*a0))*sin(phi)*	&
        sp_dbesselJ_beta0_a0(n)*Pn1(n,ct)/st-	&
        an* ((1.0,0.0)/(w*eps0*z0*a0))*sin(phi)	&
        *sp_besselJ_beta0_a0(n)*dPn1(n,ct)*st
     
! given the expansion coefficients, calculate the scattered E field in the far field    
           Ethetaf=Ethetaf+j*(exp(-j*beta0_a2)/beta0)*cos(phi)*(j**n)	&
                *(  bn*st*dPn1(n,ct)	&
                   -cn*Pn1(n,ct)/st    )
     
           Ephif=Ephif+j*(exp(-j*beta0_a2)/beta0)*sin(phi)*(j**n)	&
                *(  bn*Pn1(n,ct)/st	&
                   -cn*st*dPn1(n,ct)   )
     

! Convergence testing
         oterm=sqrt(abs(bn*bn)+abs(cn*cn))
         iterm=sqrt(abs(dn*dn)+abs(fn*fn))
	 osum=osum+oterm
	 isum=isum+iterm
	 
	 if ((iterm+oterm)/(isum+osum).lt.cvg) then
           print*,'Terminated series, n=',n
	   goto 1234
         end if
	 
       end do
       
1234   continue       
       
       sigma3d=   4.0*pi*( abs(Ethetaf) *abs(Ethetaf) + abs(Ephif)   *abs(Ephif)    )

! TOTAL FIELD         
       Etheta=Etheta
       Ephi=Ephi
! SCATTERED FIELD       
!       Etheta=Etheta-Ethetai
!       Ephi=Ephi-Ephii
              
       return       
       end
