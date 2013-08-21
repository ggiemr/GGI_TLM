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
! Evaluation of Legendre polynomials. A very crude implementation which should really only be
! used for low orders. More accurate numerical methods exist and should be implemented in future.

!
! ______________________________________________
!
       function Pn0(n,x)

       real*8 Pn0
       integer n
       real*8 x
       
       integer*4 mmax,m
       integer*4 sign
       real*8 factorial,coeff1
       
       real*8 result
       integer n4
       
       n4=n
       
       mmax=int(n/2)
       
       result=0d0
       sign=1
       do m=0,mmax
         result=result+sign*(coeff1(n4-m)/(2d0**m))*(x**(n4-2*m))/  &   
     	     factorial(m)/factorial(n4-2*m)
	 sign=sign*(-1)       
       end do
       
       Pn0=result
       
       return
    
       end
!
! ______________________________________________
!
       function Pn1(n,x)

       real*8 Pn1
       integer n
       real*8 x
       
       integer*4 mmax,m
       integer sign
       real*8 factorial,coeff1
       
       real*8 result
       integer n4
       
       n4=n
       
       mmax=int(n/2)
       
       result=0d0
       sign=1
       do m=0,mmax
         if((n-2*m-1).ge.0) then
           result=result+sign*(coeff1(n4-m)/(2d0**m))*(n4-2*m)*	&
          	 (x**(n4-2*m-1))/factorial(m)/factorial(n4-2*m)
         end if
	 sign=sign*(-1)       
       end do
       
       Pn1=-sqrt(1d0-x*x)*result
       
       return
    
       end
!
! ______________________________________________
!
       function dPn0(n,x)

       real*8 dPn0
       integer n
       real*8 x
       
       integer*4 mmax,m
       integer sign
       real*8 factorial,coeff1
       
       real*8 result
       integer n4
       
       n4=n
       
       mmax=int(n/2)
       
       result=0d0
       sign=1
       do m=0,mmax
       
         if((n4-2*m-1).ge.0) then
           result=result+sign*(coeff1(n4-m)/(2d0**m))*(n4-2*m)	&
         *(x**(n4-2*m-1))/factorial(m)/factorial(n4-2*m)
         end if
	 sign=sign*(-1)       
       end do
       
       dPn0=result
       
       return
    
       end
!
! ______________________________________________
!
       function dPn1(n,x)

       real*8 dPn1
       integer n
       real*8 x
       
       integer*4 mmax,m
       integer sign
       real*8 factorial,coeff1
       
       real*8 result1,result2,coeff
       integer n4
       
       n4=n
       
       mmax=int(n/2)
       
       result1=0d0
       result2=0d0
       sign=1
       do m=0,mmax
       
       coeff=sign*(coeff1(n4-m)/ (2d0**m))   &
     	     /factorial(m)/factorial(n4-2*m)
     
         if((n4-2*m-1).ge.0) then
           result1=result1+coeff*(n4-2*m)*(x**(n4-2*m-1))
	 end if
	 
         if((n4-2*m-2).ge.0) then
          result2=result2+coeff*(n4-2*m)*(n4-2*m-1)*(x**(n4-2*m-2))
	 end if
	 
	 sign=sign*(-1)       
	 
       end do
       
       dPn1= result1*x/sqrt(1d0-x*x)-result2*sqrt(1d0-x*x)
       
       return
    
       end
!
! ______________________________________________
!
       function factorial(n)
       
       real*8 factorial
       integer*4 n
       
       integer*4 i
       real*8 result
       
       result=1.0
       do i=1,n
         result=result*dble(i)
       end do

       factorial=result
       
       return
       end
!
! ______________________________________________
!
       function coeff1(n)
       
       real*8 coeff1
       integer*4 n
       
       integer*4 i
       real*8 result
       
       result=1d0
       do i=2,n
         result=result*(dble(i)*2d0-1d0)
       end do

       coeff1=result
       
       return
       end

