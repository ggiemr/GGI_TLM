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
!+
! NAME
!     polynomial_functions
!
! DESCRIPTION
!     definitions of functions operating on polynomial types
!
! HISTORY
!
!     started 23/01/09 CJS
!
!

MODULE polynomial_functions

USE polynomial_types

CONTAINS

!
! NAME
!     allocate_polynomial
!
! DESCRIPTION
!     set up a polynomial type with a given order
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION allocate_polynomial(order) RESULT(res)

USE polynomial_types

IMPLICIT NONE

! argument types
  integer 	:: order
! Result type
  type(polynomial) :: ans
  type(polynomial) :: res

! function definition

  ans%order=order
  allocate (ans%coeff(0:ans%order))
  ans%coeff(0:ans%order)=0d0
  
  res=ans
  
END FUNCTION allocate_polynomial

!
! NAME
!     allocate_complex_polynomial
!
! DESCRIPTION
!     set up a complex polynomial type with a given order
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION allocate_complex_polynomial(order) RESULT(res)

USE polynomial_types

IMPLICIT NONE

! argument types
  integer 	:: order
! Result type
  type(complex_polynomial) :: ans
  type(complex_polynomial) :: res

! function definition

  ans%order=order
  allocate (ans%coeff(0:ans%order))
  ans%coeff(0:ans%order)=(0d0,0d0)
  
  res=ans
  
END FUNCTION allocate_complex_polynomial

!
! NAME
!     evaluate_polynomial
!
! DESCRIPTION
!     evaluate a polynomial with real coefficients and complex argument
!
! HISTORY
!
!     started 26/01/09 CJS
!

  FUNCTION evaluate_polynomial(a,s) RESULT(res)

  USE polynomial_types
  
  type(polynomial), intent(in) :: a
  complex*16     , intent(in) :: s
  
  complex*16   		    :: res

! local variables  
  complex*16 sn
  integer n

!START

  sn=(1d0,0d0)
  res=(0d0,0d0)
  do n=0,a%order
    res=res+a%coeff(n)*sn
    sn=sn*s
  end do
 
  return
  end function evaluate_polynomial
!
! NAME
!     evaluate_complex_polynomial
!
! DESCRIPTION
!     evaluate a complex polynomial with complex coefficients and complex argument
!
! HISTORY
!
!     started 26/01/09 CJS
!

  FUNCTION evaluate_complex_polynomial(a,s) RESULT(res)

  USE polynomial_types
  
  type(complex_polynomial), intent(in) :: a
  complex*16     , intent(in) :: s
  
  complex*16   		    :: res

! local variables  
  complex*16 sn
  integer n

!START

  sn=(1d0,0d0)
  res=(0d0,0d0)
  do n=0,a%order
    res=res+a%coeff(n)*sn
    sn=sn*s
  end do
 
  return
  end function evaluate_complex_polynomial

END MODULE polynomial_functions
