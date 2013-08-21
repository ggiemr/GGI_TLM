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
! Name
!     Operators -  overloading operators for polynomials
!
! Description
!     This module defines arithmetic operations for polynomials and complex polynomials
!
! Public
!     The following are defined in this module
!
!     +
!     -
!     *
!!     /
!     =
!
! History
!
!     started 23/01/09 CJS
!
MODULE polynomial_operators

IMPLICIT NONE

PRIVATE  ! everything private unless stated explicitly

PUBLIC :: OPERATOR(+), OPERATOR(-), OPERATOR(*), ASSIGNMENT(=)

! OVERLOAD + to operate on value and polynomial data type

INTERFACE OPERATOR(+)
  MODULE PROCEDURE Addpoly, Addpoly_complex
END INTERFACE

! OVERLOAD - to operate on value and polynomial data type

INTERFACE OPERATOR(-)
  MODULE PROCEDURE Subtractpoly, Subtractpoly_complex
END INTERFACE

! OVERLOAD * to operate on value and polynomial data type

INTERFACE OPERATOR(*)
  MODULE PROCEDURE Multiplypoly, Multiplypoly_complex
END INTERFACE

! OVERLOAD = to assign a real value to polynomial types

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE Assignpoly, Assignpoly2, Assignpoly_complex, Assignpoly_complex2
END INTERFACE


! ---------------------------------------------------------
! ---------------------------------------------------------

CONTAINS

!+
! NAME
!     Addpoly
!
! DESCRIPTION
!     Add polynomial types
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION Addpoly(p1,p2) RESULT(res)

USE polynomial_types
USE polynomial_functions

IMPLICIT NONE

! argument types
type(polynomial), intent(in) :: p1,p2
! Result type
type(polynomial) :: res

! local variable
type(polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)+p2%coeff(i)
end do

res=answer

END FUNCTION Addpoly

!+
! NAME
!     Addpoly_complex
!
! DESCRIPTION
!     Add complex polynomial types
!
! HISTORY
!
!     started 4/01/13 CJS
!

FUNCTION Addpoly_complex(p1,p2) RESULT(res)

USE polynomial_types
USE polynomial_functions

IMPLICIT NONE

! argument types
type(complex_polynomial), intent(in) :: p1,p2
! Result type
type(complex_polynomial) :: res

! local variable
type(complex_polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_complex_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)+p2%coeff(i)
end do

res=answer

END FUNCTION Addpoly_complex

! ---------------------------------------------------------


!+
! NAME
!     Subtractpoly
!
! DESCRIPTION
!     Subtract polynomial types
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION Subtractpoly(p1,p2) RESULT(res)

USE polynomial_types
USE polynomial_functions

IMPLICIT NONE

! argument types
type(polynomial), intent(in) :: p1,p2
! Result type
type(polynomial) :: res

! local variable
type(polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)-p2%coeff(i)
end do

res=answer

END FUNCTION Subtractpoly


!+
! NAME
!     Subtractpoly_complex
!
! DESCRIPTION
!     Subtract complex polynomial types
!
! HISTORY
!
!     started 4/1/2013 CJS
!

FUNCTION Subtractpoly_complex(p1,p2) RESULT(res)

USE polynomial_types
USE polynomial_functions

IMPLICIT NONE

! argument types
type(complex_polynomial), intent(in) :: p1,p2
! Result type
type(complex_polynomial) :: res

! local variable
type(complex_polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_complex_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)-p2%coeff(i)
end do

res=answer

END FUNCTION Subtractpoly_complex


! ---------------------------------------------------------


!+
! NAME
!     Multiplypoly
!
! DESCRIPTION
!     Multiply polynomial types
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION Multiplypoly(p1,p2) RESULT(res)

USE polynomial_types
USE polynomial_functions

IMPLICIT NONE

! argument types
type(polynomial), intent(in) :: p1,p2
! Result type
type(polynomial) :: res

! local variables
type(polynomial) :: answer

integer i1,j1,order

! function definition

order=p1%order+p2%order
answer=allocate_polynomial(order)
answer%order=order

  do i1=0,p1%order
    do j1=0,p2%order
! product of p1(i1) p2(j1)
      answer%coeff(i1+j1)=answer%coeff(i1+j1)+p1%coeff(i1)*p2%coeff(j1)
    end do
  end do  

res=answer

END FUNCTION Multiplypoly
!+
! NAME
!     Multiplypoly_complex
!
! DESCRIPTION
!     Multiply complex polynomial types
!
! HISTORY
!
!     started 4/1/2013 CJS
!

FUNCTION Multiplypoly_complex(p1,p2) RESULT(res)

USE polynomial_types
USE polynomial_functions

IMPLICIT NONE

! argument types
type(complex_polynomial), intent(in) :: p1,p2
! Result type
type(complex_polynomial) :: res

! local variables
type(complex_polynomial) :: answer

integer i1,j1,order

! function definition

order=p1%order+p2%order
answer=allocate_complex_polynomial(order)
answer%order=order

  do i1=0,p1%order
    do j1=0,p2%order
! product of p1(i1) p2(j1)
      answer%coeff(i1+j1)=answer%coeff(i1+j1)+p1%coeff(i1)*p2%coeff(j1)
    end do
  end do  

res=answer

END FUNCTION Multiplypoly_complex
 
! ---------------------------------------------------------
! NAME
!     assignpoly
!
! DESCRIPTION
!     assign polynomial types: set to real
!
! COMMENTS
!     
! HISTORY
!
!     started 23/01/09 CJS
!

SUBROUTINE Assignpoly(res,b)

USE polynomial_types

IMPLICIT NONE

! Result type
type(polynomial), intent(out):: res
! argument types
real*8	      , intent(in) :: b

integer order

! SUBROUTINE definition

  order=0
  if (.NOT. ALLOCATED (res%coeff) ) then
! this filter not yet allocated so allocate here
    res%order=order
    ALLOCATE (res%coeff(0:order))
  else if (res%order.ne.order) then
! this filter has already been allocated but to the wrong order so reallocate  
    res%order=order
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:order))
  end if
  
  res%order=0

  res%coeff(0)=b

END SUBROUTINE Assignpoly

! ---------------------------------------------------------
! NAME
!     assignpoly2
!
! DESCRIPTION
!     assign polynomial types: set to existing polynomial
!
! COMMENTS
!     
! HISTORY
!
!     started 23/01/09 CJS
!

SUBROUTINE Assignpoly2(res,a)

USE polynomial_types

IMPLICIT NONE

! Result type
type(polynomial), intent(out):: res
! argument types
type(polynomial), intent(in) :: a

! SUBROUTINE definition

  if (.NOT. ALLOCATED (res%coeff)) then
! this filter not yet allocated so allocate here
    ALLOCATE (res%coeff(0:a%order))
  else if (res%order.ne.a%order) then
! this filter has already been allocated but to the wrong order so reallocate  
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:a%order))
  end if
  
  res%order=a%order
  res%coeff(0:res%order)=a%coeff(0:a%order)

END SUBROUTINE Assignpoly2

! ---------------------------------------------------------
! NAME
!     assignpoly_complex
!
! DESCRIPTION
!     assign complex polynomial types: set to real
!
! COMMENTS
!     
! HISTORY
!
!     started 4/1/2013 CJS
!

SUBROUTINE Assignpoly_complex(res,b)

USE polynomial_types

IMPLICIT NONE

! Result type
type(complex_polynomial), intent(out):: res
! argument types
complex*16	      , intent(in) :: b

integer order

! SUBROUTINE definition

  order=0
  if (.NOT. ALLOCATED (res%coeff) ) then
! this filter not yet allocated so allocate here
    res%order=order
    ALLOCATE (res%coeff(0:order))
  else if (res%order.ne.order) then
! this filter has already been allocated but to the wrong order so reallocate  
    res%order=order
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:order))
  end if
  
  res%order=0

  res%coeff(0)=b

END SUBROUTINE Assignpoly_complex

! ---------------------------------------------------------
! NAME
!     assignpoly_complex2
!
! DESCRIPTION
!     assign complex polynomial types: set to existing polynomial
!
! COMMENTS
!     
! HISTORY
!
!     started 4/1/2013 CJS
!

SUBROUTINE Assignpoly_complex2(res,a)

USE polynomial_types

IMPLICIT NONE

! Result type
type(complex_polynomial), intent(out):: res
! argument types
type(complex_polynomial), intent(in) :: a

! SUBROUTINE definition

  if (.NOT. ALLOCATED (res%coeff)) then
! this filter not yet allocated so allocate here
    ALLOCATE (res%coeff(0:a%order))
  else if (res%order.ne.a%order) then
! this filter has already been allocated but to the wrong order so reallocate  
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:a%order))
  end if
  
  res%order=a%order
  res%coeff(0:res%order)=a%coeff(0:a%order)

END SUBROUTINE Assignpoly_complex2


END MODULE  polynomial_operators

