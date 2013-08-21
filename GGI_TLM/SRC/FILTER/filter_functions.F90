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
!FUNCTION reduced_order_Sfilter(s1,max_freq,max_error) RESULT(res)
!FUNCTION allocate_Sfilter(aorder,border) RESULT(res)
!FUNCTION allocate_Zfilter(aorder,border) RESULT(res)
!FUNCTION allocate_Zfilter_response(aorder,border) RESULT(res)
!FUNCTION s_to_z(a,T) RESULT(res)
!FUNCTION s_to_z_warp(a,T,warp_flag,warp_scale) RESULT(res)
!FUNCTION z_to_s(a) RESULT(res)
!FUNCTION evaluate_Sfilter_frequency_response(s1,f) RESULT(res)
!FUNCTION evaluate_Sfilter_PR_frequency_response(s1,f) RESULT(res)
!FUNCTION evaluate_Sfilter_PZ_frequency_response(s1,f) RESULT(res)
!FUNCTION evaluate_Zfilter_frequency_response(z1,f) RESULT(res)
!FUNCTION Convert_filter_S_to_S_PR
!FUNCTION Convert_filter_S_PR_to_S
!FUNCTION Convert_filter_S_PZ_to_S
!FUNCTION Convert_filter_S_to_S_PZ
!
! NAME
!     filter_functions
!
! DESCRIPTION
!     definitions of functions operating on filter types
!
! HISTORY
!
!     started 22/01/09 CJS
!
!

MODULE filter_functions

USE filter_types

CONTAINS

!
! NAME
!     reduced_order_Sfilter
!
! DESCRIPTION
!     If possible calculate a reduced order filter which mimics
!     the frequency response of the input filter up to a given frequency 
!
! HISTORY
!
!     started 26/01/09 CJS
!
! COMMENTS
!     Rather simple implementation looking only at the polynimial expansion
!     of the denominator and numerator of the transfer function.
!     Does not conserve 'energy' yet in any sense i.e. has the potential
!     to introduce loss or even gain into the system. 
!     Maybe better to do something based on poles and zeros at some point 
!     as this would also take into account materials and boundary conditions
!     with poles and zeros below max_freq 
!

FUNCTION reduced_order_Sfilter(s1,max_freq,max_error) RESULT(res)

USE filter_types
USE polynomial_functions
USE constants

IMPLICIT NONE

! argument types
  type(Sfilter) :: s1
  real*8 max_freq
  real*8 max_error
  
! Result type
  type(Sfilter) :: ans
  type(Sfilter) :: res
  
! local variables
  integer den_order,num_order
  real*8  w,error
  complex*16 jw,jwn
  complex*16 den,num,num2,den2
  integer order


! START function definition

  write(*,*)'Called reduced_order_Sfilter'
  write(*,*)'max_freq=',max_freq,' max_error=',max_error

! numerator first
! calculate result at the maximum frequency of interest

  w=2d0*pi*max_freq/s1%wnorm
  jw=j*w
  num=evaluate_polynomial(s1%a,jw)
  den=evaluate_polynomial(s1%b,jw)
    
! calculate the response term by term and stop when we are closer than
! max_error. Numerator first
  write(*,*)''
  write(*,*)'Numerator'
  write(*,*)'      order         error                 max_error'

  num2=(0d0,0d0)
  jwn=(1d0,0d0)
  do order=0,s1%a%order
    num2=num2+jwn*s1%a%coeff(order)
    error=abs(num2-num)/abs(num2+num)
    write(*,*)order,error,max_error
    if (error.lt.max_error) goto 10 ! exit loop if we are close enough
    jwn=jwn*jw
  end do
  
  order=order-1
  
10 continue
  num_order=order 
  
! now set denominator - cannot be lower order than the numerator... 
  write(*,*)''
  write(*,*)'Denominator'
  write(*,*)'      order         error                 max_error'
  
  den2=(0d0,0d0)
  jwn=(1d0,0d0)
  do order=0,s1%b%order
    den2=den2+jwn*s1%b%coeff(order)
    error=abs(den2-den)/abs(den2+den)
    write(*,*)order,error,max_error
    if ((error.lt.max_error).and.(order.ge.num_order)) goto 20 ! exit loop if we are close enough
    jwn=jwn*jw
  end do
  
  order=order-1
  
20 continue
  den_order=order 

! set the values of the output filter
  ans%wnorm=s1%wnorm
  ans%a%order=num_order
  allocate (ans%a%coeff(0:ans%a%order))
  ans%b%order=den_order
  allocate (ans%b%coeff(0:ans%b%order))
  ans%a%coeff(0:ans%a%order)=s1%a%coeff(0:ans%a%order)
  ans%b%coeff(0:ans%b%order)=s1%b%coeff(0:ans%b%order)

  res=ans

END FUNCTION reduced_order_Sfilter


!
! NAME
!     allocate_Sfilter
!
! DESCRIPTION
!     set up a Sfilter type with a given order
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION allocate_Sfilter(aorder,border) RESULT(res)

USE filter_types

IMPLICIT NONE

! argument types
  integer 	:: aorder,border
  
! Result type
  type(Sfilter) :: ans
  type(Sfilter) :: res

! function definition

  ans%wnorm=1d0
  ans%a%order=aorder
  allocate (ans%a%coeff(0:ans%a%order))
  ans%b%order=border
  allocate (ans%b%coeff(0:ans%b%order))
  ans%a%coeff(0:ans%a%order)=0d0
  ans%b%coeff(0:ans%b%order)=0d0

  res=ans

END FUNCTION allocate_Sfilter

!
! NAME
!     allocate_Zfilter
!
! DESCRIPTION
!     set up a Zfilter type with a given order
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION allocate_Zfilter(aorder,border) RESULT(res)

USE filter_types

IMPLICIT NONE

! argument types
  integer 	:: aorder,border
! Result type
  type(Zfilter) :: ans
  type(Zfilter) :: res

! function definition

  ans%wnorm=1d0
  ans%T=0d0
  ans%a%order=aorder
  allocate (ans%a%coeff(0:ans%a%order))
  ans%b%order=border
  allocate (ans%b%coeff(0:ans%b%order))
  ans%a%coeff(0:ans%a%order)=0d0
  ans%b%coeff(0:ans%b%order)=0d0

  res=ans

END FUNCTION allocate_Zfilter

!
! NAME
!     allocate_Zfilter_response
!
! DESCRIPTION
!     set up a Zfilter_response type with a given order
!
! HISTORY
!
!     started 28/01/09 CJS
!

FUNCTION allocate_Zfilter_response(aorder,border) RESULT(res)

USE filter_types

IMPLICIT NONE

! argument types
  integer 	:: aorder,border
! Result type
  type(Zfilter_response) :: res

! local
  integer max_order

! function definition

  max_order=aorder
  if (border.gt.max_order) max_order=border

  res%f=0d0
  res%order=max_order
  allocate (res%w(0:res%order))
  res%w(0:res%order)=0d0


END FUNCTION allocate_Zfilter_response

!
! NAME
!     s_to_z
!
! DESCRIPTION
!     s_to_z transformation
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION s_to_z(a,T) RESULT(res)

USE filter_types

IMPLICIT NONE

! argument types
type(Sfilter), intent(in) :: a
real*8       , intent(in) :: T

! Result type
type(Zfilter) :: ans
type(Zfilter) :: res

real*8			  :: norm
integer			  :: max_order

! function definition

  max_order=a%a%order
  if (a%b%order.gt.max_order) max_order=a%b%order

  ans=allocate_Zfilter(max_order,max_order)
  call bilinear_s_to_z(a,T,ans)
  res=ans
  
  RETURN

END FUNCTION s_to_z


!
! NAME
!     s_to_z
!
! DESCRIPTION
!     s_to_z transformation
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION s_to_z_warp(a,T,warp_flag,warp_scale) RESULT(res)

USE filter_types

IMPLICIT NONE

! argument types
type(Sfilter), intent(in) :: a
real*8       , intent(in) :: T
logical, intent(in)  :: warp_flag
real*8		:: warp_scale

! Result type
type(Zfilter) :: ans
type(Zfilter) :: res

real*8			  :: norm
real*8			  :: T_scale
integer			  :: max_order

! function definition

  T_scale=T/warp_scale

  if (warp_flag) then ! warp flag set so use bicubic s_to_z

    max_order=a%a%order
    if (a%b%order.gt.max_order) max_order=a%b%order
    max_order=3*max_order
    ans=allocate_Zfilter(max_order,max_order)
    CALL bicubic_s_to_z(a,T_scale,ans)
    res=ans
  
  else ! not warp flag so use bilinear s_to_z

    max_order=a%a%order
    if (a%b%order.gt.max_order) max_order=a%b%order

    ans=allocate_Zfilter(max_order,max_order)
    CALL bilinear_s_to_z(a,T_scale,ans)
    res=ans
    
  end if
  
  return

END FUNCTION s_to_z_warp

!
! NAME
!     z_to_s
!
! DESCRIPTION
!     z_to_s transformation
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION z_to_s(a) RESULT(res)

USE filter_types

IMPLICIT NONE

! argument types
type(Zfilter), intent(in) :: a

! Result type
type(Sfilter) :: ans
type(Sfilter) :: res

real*8			  :: norm
integer			  :: max_order

! function definition

  max_order=a%a%order
  if (a%b%order.gt.max_order) max_order=a%b%order

  ans=allocate_Sfilter(ans%a%order,ans%b%order)
  call bilinear_z_to_s(a,ans)
  res=ans
  
  return

  if (a%a%order.gt.1) then
    print*,'Cant cope with more than order 1 at the moment...'
    stop
  end if
  
  ans%a%order=a%a%order
  ans%b%order=a%b%order
  ans=allocate_Sfilter(ans%a%order,ans%b%order)
  
  ans%wnorm=a%wnorm
  ans%a%coeff(0)=  a%a%coeff(0)+a%a%coeff(1)
  ans%a%coeff(1)=( a%a%coeff(0)-a%a%coeff(1))*a%T/2d0

  ans%b%coeff(0)=  a%b%coeff(0)+a%b%coeff(1)
  ans%b%coeff(1)=( a%b%coeff(0)-a%b%coeff(1))*a%T/2d0
  
  norm=ans%b%coeff(0)
  if (norm.ne.0d0) then
    ans%a%coeff(:)=ans%a%coeff(:)/norm
    ans%b%coeff(:)=ans%b%coeff(:)/norm
  end if

  res=ans

END FUNCTION z_to_s

!
! NAME
!      evaluate_Sfilter_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Laplace domain filter at
! 	a specified frequency. Filter format is Rational Function
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  FUNCTION evaluate_Sfilter_frequency_response(s1,f) RESULT(res)

USE filter_types
USE polynomial_functions
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter), intent(in) :: s1
  real*8       , intent(in) :: f
  
  complex*16   		    :: res

! local variables  
  real*8 w
  complex*16 jw,jwn,num,den
  integer n

!START

    w=2d0*pi*f/s1%wnorm
    jw=j*w
    
    num=evaluate_polynomial(s1%a,jw)
    den=evaluate_polynomial(s1%b,jw)
    
    res=num/den
  
  RETURN
  END FUNCTION evaluate_Sfilter_frequency_response
!
! NAME
!      evaluate_Sfilter_PR_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Laplace domain filter at
! 	a specified frequency. Filter format is Pole-Zero
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  FUNCTION evaluate_Sfilter_PR_frequency_response(s1,f) RESULT(res)

USE filter_types
USE polynomial_functions
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter_PR), intent(in) :: s1
  real*8       , intent(in) :: f
  
  complex*16   		    :: res

! local variables  
  real*8 w
  complex*16 jw,value
  integer n

!START

    w=2d0*pi*f/s1%wnorm
    jw=j*w
    
    value=s1%C
    
    do n=1,s1%order
      
      value=value+s1%residues(n)/(jw-s1%poles(n))
      
    end do ! next term in filter function
    
    res=value
  
  RETURN
  END FUNCTION evaluate_Sfilter_PR_frequency_response
!
! NAME
!      evaluate_Sfilter_PZ_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Laplace domain filter at
! 	a specified frequency. Filter format is Pole-Zero
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  FUNCTION evaluate_Sfilter_PZ_frequency_response(s1,f) RESULT(res)

USE filter_types
USE polynomial_functions
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter_PZ), intent(in) 	:: s1
  real*8       , intent(in) 	:: f
  
  complex*16   		    	:: res

! local variables  
  real*8 w
  complex*16 jw,value
  integer n

!START

    w=2d0*pi*f/s1%wnorm
    jw=j*w
    
    value=cmplx(s1%G)

    do n=1,s1%order
      
      value=value*(jw-s1%zeros(n))/(jw-s1%poles(n))
      
    end do ! next term in filter function
        
    res=value
   
  RETURN
  END FUNCTION evaluate_Sfilter_PZ_frequency_response

!
! NAME
!      evaluate_Zfilter_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Z domain filter at
! 	a specified frequency
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  FUNCTION evaluate_Zfilter_frequency_response(z1,f) RESULT(res)

USE filter_types
USE polynomial_functions
USE constants
  
  IMPLICIT NONE
  
  type(Zfilter), intent(in) :: z1
  real*8       , intent(in) :: f
  
  complex*16   		    :: res

! local variables  
  real*8 w
  complex*16 jw,z,zn,num,den
  integer n

!START

    w=2d0*pi*f/z1%wnorm
    jw=j*w
    z=exp(-jw*z1%T)  ! note assume expansion in z^(-1) here
    
    num=evaluate_polynomial(z1%a,z)
    den=evaluate_polynomial(z1%b,z)
    
    res=num/den
  
  return
  end function evaluate_Zfilter_frequency_response
!
! NAME
!      Convert_filter_S_to_S_PR
!     
! DESCRIPTION
!       Convert filter from rational function format to format is Pole-Residue format
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  FUNCTION Convert_filter_S_to_S_PR(S) RESULT(res)

USE filter_types
USE polynomial_types
USE polynomial_functions
USE polynomial_operators
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter), intent(in) 	:: S
  	
  type(Sfilter_PR)		:: res

! local variables  
  
  integer 			:: order
  complex*16,allocatable	:: local_poly_roots(:)
  complex*16,allocatable	:: real_roots(:)
  complex*16,allocatable	:: complex_roots(:)
  integer			:: nreal
  integer			:: ncomplex
  type(Sfilter_PR)		:: PR_filter
    
  integer			:: i
  integer 			:: row,col
  
  integer			:: matrix_dim
  
  complex*16,allocatable	:: M(:,:)
  complex*16,allocatable	:: MI(:,:)
  complex*16,allocatable	:: LHS(:)
  complex*16,allocatable	:: RHS(:)
  
  type(complex_polynomial)	:: P1,P_term

!START

  if (S%a%order.ne.S%b%order) then
    write(*,*)'Error in Convert_filter_S_to_S_PR'
    write(*,*)'a%order should be equal to b%order'
    STOP
  end if
  
  order=S%a%order
  PR_filter%order=order
  PR_filter%wnorm=S%wnorm
  
  PR_filter%n_complex_pole_pairs=0
  PR_filter%n_complex_poles=0
  PR_filter%n_real_poles=0
    
! special case for zero order, no further action required
  if (order.eq.0) then  
! 'Constant' term, C 
    PR_filter%C=S%a%coeff(order)/S%b%coeff(order)
    res=PR_filter
    RETURN    
  end if
  
! calculate poles 
  
  ALLOCATE ( local_poly_roots(1:order) )
  ALLOCATE ( real_roots(1:order) )
  ALLOCATE ( complex_roots(1:order) )

  CALL findroots(S%b,local_poly_roots,order)
  CALL rootsort(order,local_poly_roots,real_roots,		       &
  		complex_roots,nreal,ncomplex,order)
   		
  PR_filter%n_complex_pole_pairs=ncomplex
  PR_filter%n_complex_poles=ncomplex*2
  PR_filter%n_real_poles=nreal
  
  if (ALLOCATED( PR_filter%poles )) DEALLOCATE(PR_filter%poles)
  ALLOCATE( PR_filter%poles(1:order) )
  PR_filter%poles(1:nreal)=real_roots(1:nreal)
  PR_filter%poles(nreal+1:nreal+ncomplex*2)=complex_roots(1:ncomplex*2)
  
  if (ALLOCATED( PR_filter%complex_pole )) DEALLOCATE(PR_filter%complex_pole)
  ALLOCATE( PR_filter%complex_pole(1:order) )
  PR_filter%complex_pole(1:nreal)=.FALSE.
  PR_filter%complex_pole(nreal+1:nreal+ncomplex*2)=.TRUE.
  
  DEALLOCATE ( local_poly_roots )
  DEALLOCATE ( real_roots )
  DEALLOCATE ( complex_roots )

! Allocate the matrices M  and MI plus LHS and RHS for matrix equation
  matrix_dim=order+1
  ALLOCATE( M(matrix_dim,matrix_dim) )
  ALLOCATE( MI(matrix_dim,matrix_dim) )
  ALLOCATE( LHS(matrix_dim) )
  ALLOCATE( RHS(matrix_dim) )
  
  M(1:matrix_dim,1:matrix_dim)=(0d0,0d0)
  MI(1:matrix_dim,1:matrix_dim)=(0d0,0d0)
  LHS(1:matrix_dim)=(0d0,0d0)
  RHS(1:matrix_dim)=(0d0,0d0)
  
! get the LHS as the numerator polynomial of the rational function filter 

  LHS(1:matrix_dim)=S%a%coeff(0:order)/S%b%coeff(order)

! set up P_term as a first order polynomial (used in the form -pole+x )
  P_term=allocate_complex_polynomial(1)
  P_term%coeff(1)=(1d0,0d0)
  
  do col=1,order
  
    P1=(1d0,0d0)   ! initialise polynomial
    
    do i=1,order
      if (i.ne.col) then ! don't include the pole associated with this residue
        P_term%coeff(0)=-PR_filter%poles(i)
        P1=P1*P_term
      end if
    end do
    
! put the coefficients of the polynomial into the matrix M
    do row=1,order
      M(row,col)=P1%coeff(row-1)
    end do
    
  end do ! next row of the matrix M
  
! row=(order+1), polynomial multiplying Constant term
  col=order+1
  P1=(1d0,0d0)    ! initialise polynomial
  do i=1,order
    P_term%coeff(0)=-PR_filter%poles(i)
    P1=P1*P_term
  end do
! put the coefficients of the polynomial into the matrix M
  do row=1,order+1
    M(row,col)=P1%coeff(row-1)
  end do

! Invert the matrix M and solve for the coefficients of the pole/ residue expansion
  
!  CALL csvd_invert(M,matrix_dim,matrix_dim,MI,matrix_dim)
!  CALL csvd_invert_LAPACK(M,matrix_dim,matrix_dim,MI,matrix_dim)
  CALL cinvert_Gauss_Jordan(M,matrix_dim,MI,matrix_dim)
  
  CALL cmatvmul(MI,matrix_dim,matrix_dim,LHS,matrix_dim,RHS,matrix_dim)
  
  if (ALLOCATED( PR_filter%residues )) DEALLOCATE(PR_filter%residues)
  ALLOCATE( PR_filter%residues(1:order) )
  PR_filter%residues(1:order)=RHS(1:order)
  PR_filter%C=RHS(order+1)

  DEALLOCATE( M )
  DEALLOCATE( MI )
  DEALLOCATE( LHS )
  DEALLOCATE( RHS )
 
  res=PR_filter
  
  RETURN
  END FUNCTION Convert_filter_S_to_S_PR

!
! NAME
!      Convert_filter_S_PR_to_S
!     
! DESCRIPTION
!       Convert filter from Pole-Residue format to rational function format 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  FUNCTION Convert_filter_S_PR_to_S(S_PR) RESULT(res)

USE filter_types
USE polynomial_types
USE polynomial_functions
USE polynomial_operators
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter_PR), intent(in) 	:: S_PR
  	
  type(Sfilter)		:: res

! local variables  
  
  integer 			:: order

  type(Sfilter)			:: S_filter
  
  integer			:: i,term
  
  type(complex_polynomial)	:: P1,P2,P_term
  
  real*8			:: norm

!START
  
  order=S_PR%order
  S_filter=allocate_Sfilter(order,order)
  S_filter%a%order=order
  S_filter%b%order=order
  S_filter%wnorm=S_PR%wnorm
    
! special case for zero order, no further action required
  if (order.eq.0) then  
    S_filter%a%coeff(0)=S_PR%C
    S_filter%b%coeff(0)=1d0
    res=S_filter
    RETURN    
  end if
! set up P_term as a first order polynomial (used in the form -pole+x )
  P_term=allocate_complex_polynomial(1)
  P_term%coeff(1)=(1d0,0d0)

! numerator function
  P2=(0d0,0d0)   ! initialise polynomial
  
  do term=1,order
  
    P1=S_PR%residues(term)   ! initialise polynomial
    
    do i=1,order
      if (i.ne.term) then ! don't include the pole associated with this residue
        P_term%coeff(0)=-S_PR%poles(i)
        P1=P1*P_term
      end if
    end do
    
    P2=P2+P1
    
  end do ! next residue term

! constant term  
  P1=dcmplx(S_PR%C)   ! initialise polynomial
  do i=1,order
    P_term%coeff(0)=-S_PR%poles(i)
    P1=P1*P_term
  end do
  
  P2=P2+P1
  
! put the coefficients of the polynomial into the numerator of S_filter
  S_filter%a%coeff(0:order)=P2%coeff(0:order)

! denominator function
    
  P1=(1d0,0d0)    ! initialise polynomial
  do i=1,order
    P_term%coeff(0)=-S_PR%poles(i)
    P1=P1*P_term
  end do
  
! put the coefficients of the polynomial into the denominator of S_filter
  S_filter%b%coeff(0:order)=P1%coeff(0:order)

! normalise such that S_filter%b%coeff(0)=1d0
  norm=S_filter%b%coeff(0)
  if (norm.ne.0d0) then
    S_filter%a%coeff(0:order)=S_filter%a%coeff(0:order)/norm
    S_filter%b%coeff(0:order)=S_filter%b%coeff(0:order)/norm
  else
    write(*,*)'Error in Convert_filter_S_PZ_to_S'
    write(*,*)'b(0)=0d0'
    STOP
  end if

  res=S_filter
  
  RETURN
  END FUNCTION Convert_filter_S_PR_to_S
!
! NAME
!      Convert_filter_S_PZ_to_S
!     
! DESCRIPTION
!       Convert filter from Pole-Zero format to rational function format 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  FUNCTION Convert_filter_S_PZ_to_S(S_PZ) RESULT(res)

USE filter_types
USE polynomial_types
USE polynomial_functions
USE polynomial_operators
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter_PZ), intent(in) 	:: S_PZ
  	
  type(Sfilter)		:: res

! local variables  
  
  integer 			:: order

  type(Sfilter)			:: S_filter
  
  integer			:: i
  type(complex_polynomial)	:: P1,P_term
  
  real*8			:: norm

!START
  
  order=S_PZ%order
  S_filter=allocate_Sfilter(order,order)
  S_filter%a%order=order
  S_filter%b%order=order
  S_filter%wnorm=S_PZ%wnorm
    
! special case for zero order, no further action required
  if (order.eq.0) then  
    S_filter%a%coeff(0)=S_PZ%G
    S_filter%b%coeff(0)=1d0
    res=S_filter
    RETURN    
  end if

! numerator function

! set up P_term as a first order polynomial (used in the form -pole+x )
  P_term=allocate_complex_polynomial(1)
  P_term%coeff(1)=(1d0,0d0)
    
  P1=(1d0,0d0)    ! initialise polynomial
  do i=1,order
    P_term%coeff(0)=-S_PZ%zeros(i)
    P1=P1*P_term
  end do
  
! put the coefficients multiplied by Gain of the polynomial into the numerator of S_filter
  S_filter%a%coeff(0:order)=P1%coeff(0:order)*S_PZ%G

! denominator function
    
  P1=(1d0,0d0)    ! initialise polynomial
  do i=1,order
    P_term%coeff(0)=-S_PZ%poles(i)
    P1=P1*P_term
  end do
  
! put the coefficients of the polynomial into the denominator of S_filter
  S_filter%b%coeff(0:order)=P1%coeff(0:order)

! normalise such that S_filter%b%coeff(0)=1d0
  norm=S_filter%b%coeff(0)
  if (norm.ne.0d0) then
    S_filter%a%coeff(0:order)=S_filter%a%coeff(0:order)/norm
    S_filter%b%coeff(0:order)=S_filter%b%coeff(0:order)/norm
  else
    write(*,*)'Error in Convert_filter_S_PZ_to_S'
    write(*,*)'b(0)=0d0'
    STOP
  end if
  
  res=S_filter
  
  RETURN
  END FUNCTION Convert_filter_S_PZ_to_S

!
! NAME
!      Convert_filter_S_to_S_PZ
!     
! DESCRIPTION
!       Convert filter from rational function format to format is Pole-Zero format
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  FUNCTION Convert_filter_S_to_S_PZ(S) RESULT(res)

USE filter_types
USE polynomial_functions
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter), intent(in) 	:: S
  
  type(Sfilter_PZ) 		:: res

! local variables  
  
  integer order
                                   
  complex*16,allocatable	:: local_poly_roots(:)
  complex*16,allocatable	:: real_roots(:)
  complex*16,allocatable	:: complex_roots(:)
  integer			:: nreal
  integer			:: ncomplex
  
  type(Sfilter_PZ)		:: PZ_filter
       
! START

  if (S%a%order.ne.S%b%order) then
    write(*,*)'Error in Convert_filter_S_to_S_PZ'
    write(*,*)'a%order should be equal to b%order'
    STOP
  end if
  
  order=S%a%order
  PZ_filter%order=order
  PZ_filter%wnorm=S%wnorm
  
  PZ_filter%n_complex_pole_pairs=0
  PZ_filter%n_complex_poles=0
  PZ_filter%n_real_poles=0
  
  PZ_filter%n_complex_zero_pairs=0
  PZ_filter%n_complex_zeros=0
  PZ_filter%n_real_zeros=0

! 'Gain' term, G  
  PZ_filter%G=S%a%coeff(order)/S%b%coeff(order)
  
! special case for zero order, no further action required
  if (order.eq.0) then  
    res=PZ_filter
    RETURN    
  end if
  
  ALLOCATE ( local_poly_roots(1:order) )
  ALLOCATE ( real_roots(1:order) )
  ALLOCATE ( complex_roots(1:order) )

! calculate zeros

  CALL findroots(S%a,local_poly_roots,order)
  CALL rootsort(order,local_poly_roots,real_roots,		       &
  		complex_roots,nreal,ncomplex,order)
   		
  PZ_filter%n_complex_zero_pairs=ncomplex
  PZ_filter%n_complex_zeros=ncomplex*2
  PZ_filter%n_real_zeros=nreal
  
  if (ALLOCATED( PZ_filter%zeros )) DEALLOCATE(PZ_filter%zeros)
  ALLOCATE( PZ_filter%zeros(1:order) )
  PZ_filter%zeros(1:nreal)=real_roots(1:nreal)
  PZ_filter%zeros(nreal+1:nreal+ncomplex*2)=complex_roots(1:ncomplex*2)
  
  if (ALLOCATED( PZ_filter%complex_zero )) DEALLOCATE(PZ_filter%complex_zero)
  ALLOCATE( PZ_filter%complex_zero(1:order) )
  PZ_filter%complex_zero(1:nreal)=.FALSE.
  PZ_filter%complex_zero(nreal+1:nreal+ncomplex*2)=.TRUE.
  
! calculate poles 

  CALL findroots(S%b,local_poly_roots,order)
  CALL rootsort(order,local_poly_roots,real_roots,		       &
  		complex_roots,nreal,ncomplex,order)
   		
  PZ_filter%n_complex_pole_pairs=ncomplex
  PZ_filter%n_complex_poles=ncomplex*2
  PZ_filter%n_real_poles=nreal
  
  if (ALLOCATED( PZ_filter%poles )) DEALLOCATE(PZ_filter%poles)
  ALLOCATE( PZ_filter%poles(1:order) )
  PZ_filter%poles(1:nreal)=real_roots(1:nreal)
  PZ_filter%poles(nreal+1:nreal+ncomplex*2)=complex_roots(1:ncomplex*2)
  
  if (ALLOCATED( PZ_filter%complex_pole )) DEALLOCATE(PZ_filter%complex_pole)
  ALLOCATE( PZ_filter%complex_pole(1:order) )
  PZ_filter%complex_pole(1:nreal)=.FALSE.
  PZ_filter%complex_pole(nreal+1:nreal+ncomplex*2)=.TRUE.
  
  DEALLOCATE ( local_poly_roots )
  DEALLOCATE ( real_roots )
  DEALLOCATE ( complex_roots )
  
  res=PZ_filter
  
  RETURN
  END FUNCTION Convert_filter_S_to_S_PZ
  
END MODULE filter_functions
