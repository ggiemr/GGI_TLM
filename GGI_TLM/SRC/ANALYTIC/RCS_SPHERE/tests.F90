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
       
       subroutine test_legendre()
       
USE constants

       implicit none
       
! Legendre polynomials
       real*8 Pn0,Pn1,dPn0,dPn1
       
     
       real*8 x1,x2,x3
       real*8 f1,f2,f3
       real*8 dfdx
     
       real*8 x,dx
       
       real*8 theta
       
       integer n
       
       print*,'Test Legendre polynomial evaluation...'
       
       do n=5,6
       
         x=0.5
         dx=1e-6
       
         x1=x-dx
         x2=x
         x3=x+dx
       
         f1=Pn0(n,x1)
         f2=Pn0(n,x2)
         f3=Pn0(n,x3)
       
         dfdx=dPn0(n,x2)
       
         print*,' '
       
         print*,'n=',n
         print*,'x=',x
         print*,'Pn0(n,x)=',f1
         print*,'Pn0(n,x)=',f2
         print*,'Pn0(n,x)=',f3
         print*,'analytic  dPn0/dx=',dfdx
         print*,'numerical dPn0/dx=',(f3-f1)/(2.0*dx)
       
       
         print*,' '
       
         f1=Pn1(n,x1)
         f2=Pn1(n,x2)
         f3=Pn1(n,x3)
       
         dfdx=dPn1(n,x2)
       
         print*,'Pn1(n,x)=',f1
         print*,'Pn1(n,x)=',f2
         print*,'Pn1(n,x)=',f3
         print*,'analytic  dPn1/dx=',dfdx
         print*,'numerical dPn1/dx=',(f3-f1)/(2.0*dx)
       
       end do
       
       print*,' checks at theta=pi'
       theta=pi*1.00001
       do n=1,3
         print*,' '
         print*,'theta=',theta
         print*,'-Pn1(n,cos(theta))/sin(theta)',	&
              -Pn1(n,cos(theta))/sin(theta)
         print*,'((-1)**n)*(n*(n+1))/2.0      ',((-1)**n)*(n*(n+1))/2.0
         print*,'dPn1(n,cos(theta))*sin(theta)',	&
              dPn1(n,cos(theta))*sin(theta)
         print*,'-((-1)**n)*(n*(n+1))/2.0     ',-((-1)**n)*(n*(n+1))/2.0
         print*,' '
       end do
       
       print*,'End of Legendre Polynomial evaluation checks'
       print*,' '
       
       return       
       end
!
!_____________________________________________________
!       
       
       subroutine test_bessel()

! Note: Still need to set up large argument test
       
USE constants

       implicit none
                     
       integer,parameter :: nmax=2
       
! bessel functions    
   
       complex*16 sp_besselJ(0:nmax,1:5)
       complex*16 sp_bessely(0:nmax,1:5)
       complex*16 sp_besselH1(0:nmax,1:5)
       complex*16 sp_besselH2(0:nmax,1:5)
       complex*16 sp_dbesselJ(0:nmax,1:5)
       complex*16 sp_dbessely(0:nmax,1:5)
       complex*16 sp_dbesselH1(0:nmax,1:5)
       complex*16 sp_dbesselH2(0:nmax,1:5)
               
       complex*16 sp_besselJ_temp(0:nmax)
       complex*16 sp_bessely_temp(0:nmax)
       complex*16 sp_dbesselJ_temp(0:nmax)
       complex*16 sp_dbessely_temp(0:nmax)
     
       complex*16 x(5)
       complex*16 f1,f2,f3
       complex*16 dfdx
     
       complex*16 x0,dx
       
       integer n,i,constantJ,constantY
       
       print*,'Test Bessel evaluation...'
       
       x0=(0.5,0.2)
       dx=(1e-6,0.0)
       
       x(1)=x0-dx
       x(2)=x0
       x(3)=x0+dx
       x(4)=(1d-5,5d-6)
       x(5)=(1000D0,500D0)
       
       do i=1,5
       
         CALL spherical_bessel( x(i), nmax, sp_besselJ_temp,	&
                                            sp_dbesselJ_temp,	&
					    sp_bessely_temp,	&
					    sp_dbessely_temp )
					
         sp_besselJ(0:nmax,i)= sp_besselJ_temp(0:nmax)     
         sp_besselY(0:nmax,i)= sp_besselY_temp(0:nmax)     
         sp_dbesselJ(0:nmax,i)= sp_dbesselJ_temp(0:nmax)     
         sp_dbesselY(0:nmax,i)= sp_dbesselY_temp(0:nmax) 
	     
         sp_besselH1(0:nmax,i)= sp_besselJ_temp(0:nmax)+j*sp_besselY_temp(0:nmax)
         sp_besselH2(0:nmax,i)= sp_besselJ_temp(0:nmax)-j*sp_besselY_temp(0:nmax)

         sp_dbesselH1(0:nmax,i)= sp_dbesselJ_temp(0:nmax)+j*sp_dbesselY_temp(0:nmax)
         sp_dbesselH2(0:nmax,i)= sp_dbesselJ_temp(0:nmax)-j*sp_dbesselY_temp(0:nmax)
	  
       end do
              
       do n=0,2
       
         print*,' '       
         print*,'n=',n
         print*,'x=',x0
	 
	 
         f1=                sp_besselJ(n,1)
         f2=                sp_besselJ(n,2)
         f3=                sp_besselJ(n,3)	  
         dfdx=             sp_dbesselJ(n,2)	  
         print*,' '       
         print*,'           sp_besselJ(n,x)=',f2
         print*,'analytic  dsp_besselJ/dx=  ',dfdx
         print*,'numerical dsp_besselJ/dx=  ',(f3-f1)/(2.0*dx)
	 
         f1=                sp_besselY(n,1)
         f2=                sp_besselY(n,2)
         f3=                sp_besselY(n,3)	  
         dfdx=             sp_dbesselY(n,2)	  
         print*,' '       
         print*,'           sp_besselY(n,x)=',f2
         print*,'analytic  dsp_besselY/dx=  ',dfdx
         print*,'numerical dsp_besselY/dx=  ',(f3-f1)/(2.0*dx)
	 
         f1=                sp_besselH1(n,1)
         f2=                sp_besselH1(n,2)
         f3=                sp_besselH1(n,3)	   
         dfdx=             sp_dbesselH1(n,2)	   
         print*,' '       
         print*,'           sp_besselH1(n,x)=',f2
         print*,'analytic  dsp_besselH1/dx=  ',dfdx
         print*,'numerical dsp_besselH1/dx=  ',(f3-f1)/(2.0*dx)
	 
         f1=                sp_besselH2(n,1)
         f2=                sp_besselH2(n,2)
         f3=                sp_besselH2(n,3)	   
         dfdx=             sp_dbesselH2(n,2)	   
         print*,' '       
         print*,'           sp_besselH2(n,x)=',f2
         print*,'analytic  dsp_besselH2/dx=  ',dfdx
         print*,'numerical dsp_besselH2/dx=  ',(f3-f1)/(2.0*dx)
	 	                      
       end do

!*****************************************                     
       
       print*,' '
       print*,'Wronskian checks'
       do n=0,2
       
         print*,'n=',n
	 
         print*,sp_besselJ(n,2) *sp_dbesselY(n,2)-	&
                sp_besselY(n,2) *sp_dbesselJ(n,2)
         print*,(1d0,0d0)/(x(2)*x(2))
	 print*,' '
         print*,sp_besselH1(n,2) *sp_dbesselH2(n,2)-	&
                sp_besselH2(n,2) *sp_dbesselH1(n,2)
         print*,(0d0,-2d0)/(x(2)*x(2))
	 
       end do
!*****************************************              
       
       print*,' '
       print*,'Small argument checks'
       print*,'x   =',x(4)
       
       do n=0,2
       
         constantJ=1
	 constantY=-1
	 do i=1,n
	   constantJ=constantJ*(2*i+1)
	   constantY=constantY*(2*i-1)
	 end do
       
         print*,'n   =',n
         print*,'J(x)=',sp_besselJ(n,4)
         print*,'    =',(x(4)**n)/constantJ
	 
         print*,'Y(x)=',sp_besselY(n,4)
         print*,'    =',constantY/(x(4)**(n+1))
	 
       end do

!*****************************************                     
       
!       print*,' '
!       print*,'Large argument checks'
!       
!       print*,'x   =',x(5)
!       do n=0,2
!         print*,'n   =',n
!         print*,'H2(x)=            ',sp_besselH2(n,5)
!         print*,'     ='
!         print*,'dH2(x)=           ',sp_dbesselH2(n,5)
!         print*,'      ='
!       end do
       
       return       
       end
