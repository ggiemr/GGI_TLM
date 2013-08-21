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
! NAME
!    find_roots
!
! DESCRIPTION
!       find all the roots of a polynomial using an eigenvalue method
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 17/05/05 CJS
!     Replace Laguer's method with an eigenvalue method 30/1/2013
!
! COMMENTS
! 
! Altered to make arrays dynamic.
! Assume the array roots is already allocated
! 

subroutine findroots(a,roots,poly_order)

! Modules used

USE polynomial_types

IMPLICIT NONE

!Argument list variables

type(polynomial), intent(in) 	:: a
integer     			:: poly_order
complex*16  			:: roots(poly_order)

! Local variables

  integer*4 	:: n
  integer*4 	:: low
  integer*4 	:: igh
  real*8  	:: h(poly_order,poly_order)
  real*8  	:: wi(poly_order)
  real*8  	:: wr(poly_order)
  integer*4  	:: ierr
  
  integer	:: i,j
  real*8	:: eps,test1,test2
  
! START


! CALL EISPACK subroutine to calculate the eigenvalues of an upper Hessenberg matrix

! Fill the uppper Hessenberg matrix
  h(1:poly_order,1:poly_order)=0d0
  
  do i=1,poly_order
    h(1,i)=-a%coeff(poly_order-i)/a%coeff(poly_order)
    if (i.ne.poly_order) then
      h(i+1,i)=1d0
    end if
  end do

  n=poly_order
  low=1
  igh=n
  
  CALL hqr ( n, low, igh, h, wr, wi, ierr )
  
  if (ierr.ne.0) then
    write(*,*)'Error in findroots'
    write(*,*)'Not all roots of the polynomial could be found correctly'
    STOP
  end if
    
  eps=1d-10

  do i=1,poly_order
  
! look for very small imaginary part and remove - not sure if this is required now we use eispack...
    test1=abs( wi(i) )
    test2=abs( eps*wr(i) )
    
    if (test1.le.test2) wi(i)=0d0
    
    roots(i)=dcmplx( wr(i),wi(i) )

  end do

! END
 
return
end


!
! ______________________________________________________________
!
!
! NAME
!    root_sort
!
! DESCRIPTION
!   
!    Sort real and complex roots
!
! SEE ALSO
!  
!
! HISTORY
!
!     started 17/05/05 CJS
!
! COMMENTS
!     needs some tidying up...

 subroutine rootsort(ordert,roots,rroots,             &
                           croots,nreal,ncomplex,maxordert)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer maxordert
 complex*16 roots(maxordert)
 complex*16 rroots(maxordert),croots(maxordert)
 integer ordert,nreal,ncomplex
       
!Local variables

 complex*16 swap,sum
 real*8 tangle,tanglesum
 integer i,roottype(maxordert),fpair
 integer k
 
! START

 do i=1,ordert
   if (abs(roots(i)).eq.0d0) then 
     roottype(i)=-1
   else if (abs(dble((0d0,-1d0)*roots(i))/abs(roots(i)))    &
					  .gt.1d-8) then
!  complex root found
     roottype(i)=1
   else
! real root found
     roottype(i)=-1
   end if  

 end do
!
!  reorder roots, real roots first then complex roots in pairs
!  
 nreal=0
 ncomplex=0
 do i=1,ordert
   if (roottype(i).eq.-1) then
! real root found
     nreal=nreal+1
     rroots(nreal)=roots(i)
     roottype(i)=0
   else if (roottype(i).eq.1) then
! complex root found
     ncomplex=ncomplex+1
     croots(ncomplex)=roots(i)
     roottype(i)=0
   end if
 end do 
!
 ncomplex=ncomplex/2

! pair up complex roots

 do i=1,ncomplex
   fpair=(i-1)*2+1
   if (dble(croots(fpair)).ne.0d0) then
     tangle=abs(dble((0d0,-1d0)*croots(fpair)/        &
		   dble(croots(fpair))))
   else
     tangle=1d10
   end if
   do k=fpair+1,ncomplex*2
     sum=croots(fpair)+croots(k)
     if (abs(sum).ne.0d0) then
       tanglesum=abs(dble((0d0,-1d0)*sum/	      &
			     abs(sum)))
     else
       tanglesum=0d0
     end if
     if (tanglesum/tangle.lt.1d-5) then
! conjugate root found so shift next to first root of the pair
       swap=croots(fpair+1)
       croots(fpair+1)=croots(k)
       croots(k)=swap
     end if
   end do

 end do

 do i=1,nreal
   roots(i)=rroots(i)
 end do
 do i=1,ncomplex*2
   roots(nreal+i)=croots(i)
 end do

! END

 end
