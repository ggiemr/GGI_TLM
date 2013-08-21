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
       
       subroutine rcs_coated(theta,phi,w,a1,a2,eps2,mu2,nmax,Etheta,Ephi,Htheta,Hphi,sigma3d)

USE constants

       implicit none
              
       real*8 theta,phi,w,a1,a2
       complex*16 eps2,mu2
       integer nmax
       complex*16 Etheta,Ephi,Htheta,Hphi
       complex*16 Ethetao,Ephio,Hthetao,Hphio
       real*8 sigma3d
       real*8 r
        
! bessel functions       
       
       complex*16 sp_besselJ_beta0_r(0:nmax)
       complex*16 sp_bessely_beta0_r(0:nmax)
       complex*16 sp_besselH1_beta0_r(0:nmax)
       complex*16 sp_besselH2_beta0_r(0:nmax)
       complex*16 sp_dbesselJ_beta0_r(0:nmax)
       complex*16 sp_dbessely_beta0_r(0:nmax)
       complex*16 sp_dbesselH1_beta0_r(0:nmax)
       complex*16 sp_dbesselH2_beta0_r(0:nmax)
       
       complex*16 sp_besselJ_beta0_a2(0:nmax)
       complex*16 sp_bessely_beta0_a2(0:nmax)
       complex*16 sp_besselH1_beta0_a2(0:nmax)
       complex*16 sp_besselH2_beta0_a2(0:nmax)
       complex*16 sp_dbesselJ_beta0_a2(0:nmax)
       complex*16 sp_dbessely_beta0_a2(0:nmax)
       complex*16 sp_dbesselH1_beta0_a2(0:nmax)
       complex*16 sp_dbesselH2_beta0_a2(0:nmax)
       
       complex*16 sp_besselJ_beta2_a1(0:nmax)
       complex*16 sp_bessely_beta2_a1(0:nmax)
       complex*16 sp_besselH1_beta2_a1(0:nmax)
       complex*16 sp_besselH2_beta2_a1(0:nmax)
       complex*16 sp_dbesselJ_beta2_a1(0:nmax)
       complex*16 sp_dbessely_beta2_a1(0:nmax)
       complex*16 sp_dbesselH1_beta2_a1(0:nmax)
       complex*16 sp_dbesselH2_beta2_a1(0:nmax)
       
       complex*16 sp_besselJ_beta2_a2(0:nmax)
       complex*16 sp_bessely_beta2_a2(0:nmax)
       complex*16 sp_besselH1_beta2_a2(0:nmax)
       complex*16 sp_besselH2_beta2_a2(0:nmax)
       complex*16 sp_dbesselJ_beta2_a2(0:nmax)
       complex*16 sp_dbessely_beta2_a2(0:nmax)
       complex*16 sp_dbesselH1_beta2_a2(0:nmax)
       complex*16 sp_dbesselH2_beta2_a2(0:nmax)
       
       complex*16 sp_besselJ_beta2_r(0:nmax)
       complex*16 sp_bessely_beta2_r(0:nmax)
       complex*16 sp_besselH1_beta2_r(0:nmax)
       complex*16 sp_besselH2_beta2_r(0:nmax)
       complex*16 sp_dbesselJ_beta2_r(0:nmax)
       complex*16 sp_dbessely_beta2_r(0:nmax)
       complex*16 sp_dbesselH1_beta2_r(0:nmax)
       complex*16 sp_dbesselH2_beta2_r(0:nmax)
       
! Legendre polynomials
       real*8 Pn0,Pn1,dPn1
       
! constants       
       real*8 lambda
       
       complex*16 beta0,beta2
       complex*16 beta0_a2,beta2_a2,beta2_a1,beta2_r,beta0_r
       
       complex*16 z2
       
! matrix elements
       complex*16 k1d,k1e       
       complex*16 k2f,k2g       
       complex*16 k3b,k3d,k3e      
       complex*16 k4c,k4f,k4g      
       complex*16 k5b,k5d,k5e      
       complex*16 k6c,k6f,k6g
       complex*16 r3,r4,r5,r6
       
       complex*16 a11,a12,a13,a21,a22,a23,a31,a32,a33
       complex*16 u11,u12,u21,u22,det,rdet
       
       complex*16 an,bn,cn,dn,en,fn,gn
       
       complex*16 Ethetai,Ephii
       complex*16 Ethetaf,Ephif

       real*8 ct,st
       
       integer n
       
       real*8 oterm,iterm,osum,isum,norm,cvg
       
! START
       
       cvg=1d-6      ! convergence criterion
       
       osum=0D0
       isum=0D0
       
       r=a2

! set parameters       

       beta0=cmplx(w/c0)
       
       beta2=w*sqrt(mu2*eps2)
       if (aimag(beta2).gt.0D0) beta2=-beta2
       z2=sqrt(mu2/eps2)
       
       beta0_a2=beta0*a2
       beta2_a1=beta2*a1
       beta2_a2=beta2*a2
       beta2_r=beta2*r
       beta0_r=beta0*r
       
! evaluate bessel functions

       CALL spherical_bessel2( beta0_r, nmax, sp_besselJ_beta0_r,	&
                                             sp_dbesselJ_beta0_r,	&
					     sp_bessely_beta0_r,	&
					     sp_dbessely_beta0_r )
					   
       sp_besselH1_beta0_r(0:nmax)=sp_besselJ_beta0_r(0:nmax)+j*sp_bessely_beta0_r(0:nmax)
       sp_besselH2_beta0_r(0:nmax)=sp_besselJ_beta0_r(0:nmax)-j*sp_bessely_beta0_r(0:nmax)
       sp_dbesselH1_beta0_r(0:nmax)=sp_dbesselJ_beta0_r(0:nmax)+j*sp_dbessely_beta0_r(0:nmax)
       sp_dbesselH2_beta0_r(0:nmax)=sp_dbesselJ_beta0_r(0:nmax)-j*sp_dbessely_beta0_r(0:nmax)

       CALL spherical_bessel2( beta0_a2, nmax, sp_besselJ_beta0_a2,	&
                                             sp_dbesselJ_beta0_a2,	&
					     sp_bessely_beta0_a2,	&
					     sp_dbessely_beta0_a2 )
					   
       sp_besselH1_beta0_a2(0:nmax)=sp_besselJ_beta0_a2(0:nmax)+j*sp_bessely_beta0_a2(0:nmax)
       sp_besselH2_beta0_a2(0:nmax)=sp_besselJ_beta0_a2(0:nmax)-j*sp_bessely_beta0_a2(0:nmax)
       sp_dbesselH1_beta0_a2(0:nmax)=sp_dbesselJ_beta0_a2(0:nmax)+j*sp_dbessely_beta0_a2(0:nmax)
       sp_dbesselH2_beta0_a2(0:nmax)=sp_dbesselJ_beta0_a2(0:nmax)-j*sp_dbessely_beta0_a2(0:nmax)

       CALL spherical_bessel2( beta2_a1, nmax, sp_besselJ_beta2_a1,	&
                                             sp_dbesselJ_beta2_a1,	&
					     sp_bessely_beta2_a1,	&
					     sp_dbessely_beta2_a1 )
					   
       sp_besselH1_beta2_a1(0:nmax)=sp_besselJ_beta2_a1(0:nmax)+j*sp_bessely_beta2_a1(0:nmax)
       sp_besselH2_beta2_a1(0:nmax)=sp_besselJ_beta2_a1(0:nmax)-j*sp_bessely_beta2_a1(0:nmax)
       sp_dbesselH1_beta2_a1(0:nmax)=sp_dbesselJ_beta2_a1(0:nmax)+j*sp_dbessely_beta2_a1(0:nmax)
       sp_dbesselH2_beta2_a1(0:nmax)=sp_dbesselJ_beta2_a1(0:nmax)-j*sp_dbessely_beta2_a1(0:nmax)

       CALL spherical_bessel2( beta2_a2, nmax, sp_besselJ_beta2_a2,	&
                                             sp_dbesselJ_beta2_a2,	&
					     sp_bessely_beta2_a2,	&
					     sp_dbessely_beta2_a2 )
					   
       sp_besselH1_beta2_a2(0:nmax)=sp_besselJ_beta2_a2(0:nmax)+j*sp_bessely_beta2_a2(0:nmax)
       sp_besselH2_beta2_a2(0:nmax)=sp_besselJ_beta2_a2(0:nmax)-j*sp_bessely_beta2_a2(0:nmax)
       sp_dbesselH1_beta2_a2(0:nmax)=sp_dbesselJ_beta2_a2(0:nmax)+j*sp_dbessely_beta2_a2(0:nmax)
       sp_dbesselH2_beta2_a2(0:nmax)=sp_dbesselJ_beta2_a2(0:nmax)-j*sp_dbessely_beta2_a2(0:nmax)

       CALL spherical_bessel2( beta2_r, nmax, sp_besselJ_beta2_r,	&
                                             sp_dbesselJ_beta2_r,	&
					     sp_bessely_beta2_r,	&
					     sp_dbessely_beta2_r )
					   
       sp_besselH1_beta2_r(0:nmax)=sp_besselJ_beta2_r(0:nmax)+j*sp_bessely_beta2_r(0:nmax)
       sp_besselH2_beta2_r(0:nmax)=sp_besselJ_beta2_r(0:nmax)-j*sp_bessely_beta2_r(0:nmax)
       sp_dbesselH1_beta2_r(0:nmax)=sp_dbesselJ_beta2_r(0:nmax)+j*sp_dbessely_beta2_r(0:nmax)
       sp_dbesselH2_beta2_r(0:nmax)=sp_dbesselJ_beta2_r(0:nmax)-j*sp_dbessely_beta2_r(0:nmax)
       
       lambda=2D0*pi*c0/w

! reset fields
       Ethetai=(0d0,0d0)
       Ephii=(0d0,0d0)
       Ethetaf=(0d0,0d0)
       Ephif=(0d0,0d0)
       
       Etheta=(0d0,0d0)
       Ephi=(0d0,0d0)
       Htheta=(0d0,0d0)
       Hphi=(0d0,0d0)
       
       Ethetao=(0d0,0d0)
       Ephio=(0d0,0d0)
       Hthetao=(0d0,0d0)
       Hphio=(0d0,0d0)

       do n=1,nmax

! incident field expansion coefficient       
         an=(j**(-n))*(2.0*n+1)/(n*(n+1))
	 
! matrix elements for solution of field expansion coefficients
         k1d=sp_dbesselH1_beta2_a1(n)
         k1e=sp_dbesselH2_beta2_a1(n)
	 
         k2f=sp_besselH1_beta2_a1(n)
         k2g=sp_besselH2_beta2_a1(n)
	 
	 k3b=( beta0/(mu0*eps0))*sp_dbesselH2_beta0_a2(n)
	 k3d=(-beta2/(mu2*eps2))*sp_dbesselH1_beta2_a2(n)
	 k3e=(-beta2/(mu2*eps2))*sp_dbesselH2_beta2_a2(n)
	 
	 k4c=(-1D0/(eps0*z0))*sp_besselH2_beta0_a2(n)
	 k4f=( 1D0/(eps2*z2))*sp_besselH1_beta2_a2(n)
	 k4g=( 1D0/(eps2*z2))*sp_besselH2_beta2_a2(n)
	 
	 k5b=(-1D0/mu0)*sp_besselH2_beta0_a2(n)
	 k5d=( 1D0/mu2)*sp_besselH1_beta2_a2(n)
	 k5e=( 1D0/mu2)*sp_besselH2_beta2_a2(n)
	 
	 k6c=( beta0/(mu0*eps0*z0))*sp_dbesselH2_beta0_a2(n)
	 k6f=(-beta2/(mu2*eps2*z2))*sp_dbesselH1_beta2_a2(n)
	 k6g=(-beta2/(mu2*eps2*z2))*sp_dbesselH2_beta2_a2(n)

! RHS elements for solution of field expansion coefficients

         r3=-an*( beta0/(mu0*eps0))*sp_dbesselJ_beta0_a2(n)
	 r4=-an*(-1D0/(eps0*z0))*sp_besselJ_beta0_a2(n)
	 r5=-an*(-1D0/mu0)*sp_besselJ_beta0_a2(n)
	 r6=-an*( beta0/(mu0*eps0*z0))*sp_dbesselJ_beta0_a2(n)

! solution of submatrix equation 1
         a11=(0d0,0d0)
	 a12=k1d
	 a13=k1e
	 a21=k3b
	 a22=k3d
	 a23=k3e
	 a31=k5b
	 a32=k5d
	 a33=k5e
	
	 u11=a21
	 u12=a22-(a23*a12/a13)
	 u21=a31
	 u22=a32-(a33*a12/a13)
	 
	 norm=(abs(U11)+abs(U12)+abs(U21)+abs(U22))
	 
	 U11=U11/norm
	 U12=U12/norm
	 R3= R3 /norm
	 
	 U21=U21/norm
	 U22=U22/norm
	 R5= R5 /norm

         det=u11*u22-u12*u21
	 	 
	 if (abs(det).ne.0D0) then
	   rdet=1D0/det
	 else 
	   rdet=0D0
	 end if
	 
	 bn=( u22*r3-u12*r5)*rdet
	 dn=(-u21*r3+u11*r5)*rdet
	 en=-a12/a13*dn
	 
! solution of submatrix equation 2	 
         a11=(0d0,0d0)
	 a12=k2f
	 a13=k2g
	 a21=k4c
	 a22=k4f
	 a23=k4g
	 a31=k6c
	 a32=k6f
	 a33=k6g
	
	 u11=a21
	 u12=a22-(a23*a12/a13)
	 u21=a31
	 u22=a32-(a33*a12/a13)
	 
	 norm=(abs(U11)+abs(U12)+abs(U21)+abs(U22))
	 
	 U11=U11/norm
	 U12=U12/norm
	 R4= R4 /norm
	 
	 U21=U21/norm
	 U22=U22/norm
	 R6= R6 /norm

         det=u11*u22-u12*u21
	 	 
	 if (abs(det).ne.0D0) then
	   rdet=1D0/det
	 else 
	   rdet=0D0
	 end if
	 	 
	 cn=( u22*r4-u12*r6)*rdet
	 fn=(-u21*r4+u11*r6)*rdet
	 gn=-a12/a13*fn
	 
! cope with the limiting case of theta=0,180 degrees...     
	 
! cope with the limiting case of theta=0,180 degrees...     
	 ct=cos(theta)	 
	 if(ct.eq.1.0) ct=1.0-0.0000001
	 if(ct.eq.-1.0) ct=-1.0+0.0000001
         st=sqrt(1.0-ct*ct)

! near field expression: total field in dielectric sphere

         Etheta=Etheta	&
     	 -(beta2/(j*w*w*mu2*eps2*r))*cos(phi)*  &
     	 (dn*sp_dbesselH1_beta2_r(n)+en*sp_dbesselH2_beta2_r(n))        &
     	 *dPn1(n,ct)*st-        &
     	 ((1D0,0D0)/(w*eps2*z2*r))*cos(phi)*    &
     	 (fn*sp_besselH1_beta2_r(n)+gn*sp_besselH2_beta2_r(n))  &
     	  *Pn1(n,ct)/st
     
         Ephi=Ephi+	&
     	  (-beta2/(j*w*w*mu2*eps2*r))*sin(phi)* &
     	 (dn*sp_dbesselH1_beta2_r(n)+en*sp_dbesselH2_beta2_r(n))        &
     	 *Pn1(n,ct)/st- &
     	  ((1D0,0D0)/(w*eps2*z2*r))*sin(phi)    &
     	 *(fn*sp_besselH1_beta2_r(n)+gn*sp_besselH2_beta2_r(n)) &
     	 *dPn1(n,ct)*st


! near field expression: field on outside of cylinder

         Ethetao=Ethetao	&
     	 -(beta0/(j*w*w*mu0*eps0*a2))*cos(phi)* &
     	 (an*sp_dbesselJ_beta0_a2(n)+bn*sp_dbesselH2_beta0_a2(n) )      &
     	  *dPn1(n,ct)*st-       &
     	 ((1D0,0D0)/(w*eps0*z0*a2))*cos(phi)*   &
     	 (an*sp_besselJ_beta0_a2(n)+cn*sp_besselH2_beta0_a2(n)  )       &
     	  *Pn1(n,ct)/st     

         Ephio=Ephio+	&
     	  (-beta0/(j*w*w*mu0*eps0*a2))*sin(phi) &
     	 *(an*sp_dbesselJ_beta0_a2(n)+bn*sp_dbesselH2_beta0_a2(n)   )   &
     	 *Pn1(n,ct)/st- &
     	  ((1D0,0D0)/(w*eps0*z0*a2))*sin(phi)   &
     	 *(an*sp_besselJ_beta0_a2(n)+cn*sp_besselH2_beta0_a2(n)  )      &
     	 *dPn1(n,ct)*st     


! near field expression: incident field

         Ethetai=Ethetai	&
     	 -an* (beta0/(j*w*w*mu0*eps0*r))*cos(phi)*      &
     	  sp_dbesselJ_beta0_r(n)*dPn1(n,ct)*st- &
     	 an* ((1D0,0D0)/(w*eps0*z0*r))*cos(phi)*        &
     	 sp_besselJ_beta0_r(n)*Pn1(n,ct)/st     

         Ephii=Ephii+	&
      	 an* (-beta0/(j*w*w*mu0*eps0*r))*sin(phi)*      &
      	 sp_dbesselJ_beta0_r(n)*Pn1(n,ct)/st-   &
      	 an* ((1D0,0D0)/(w*eps0*z0*r))*sin(phi) &
      	 *sp_besselJ_beta0_r(n)*dPn1(n,ct)*st	 

! given the expansion coefficients, calculate the scattered E field in region 1     
           Ethetaf=Ethetaf+j*(exp(-j*beta0_a2)/beta0)*cos(phi)*(j**n)	&
                 *(  bn*st*dPn1(n,ct)	&
                    -cn*Pn1(n,ct)/st    )
     
           Ephif=Ephif+j*(exp(-j*beta0_a2)/beta0)*sin(phi)*(j**n)	&
                 *(  bn*Pn1(n,ct)/st	&
                    -cn*st*dPn1(n,ct)   )

! Convergence testing
         oterm=sqrt(abs(bn*bn)+abs(cn*cn))
         iterm=sqrt(abs(dn*dn)+abs(en*en)+abs(fn*fn)+abs(gn*gn))
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
!       Etheta=Ethetao-Ethetai
!       Ephi=Ephio-Ephii
              
       return       
       end
