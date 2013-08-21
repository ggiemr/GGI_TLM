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
!     filter_operators -  overloading operators for filter functions
!
! Description
!     This module defines arithmetic operations with filter functions
!
! Public
!     The following are defined in this module
!
!     *
!     /
!     =
!
! History
!
!     started 22/01/09 CJS
!
MODULE filter_operators

IMPLICIT NONE

PRIVATE  ! everything private unless stated explicitly

PUBLIC :: OPERATOR(*), OPERATOR(/),ASSIGNMENT(=)

! OVERLOAD * to operate on value and filter data type

INTERFACE OPERATOR(*)
  MODULE PROCEDURE MultiplySfilter, MultiplyZfilter
END INTERFACE

! OVERLOAD / to operate on value and filter data type

INTERFACE OPERATOR(/)
  MODULE PROCEDURE DivideSfilter, DivideZfilter
END INTERFACE

! OVERLOAD = to assign filter types

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE AssignSfilter, AssignZfilter, AssignSfilter2, 	&
                   AssignZfilter2, AssignS_PZ_filter, AssignS_PR_filter
END INTERFACE


! ---------------------------------------------------------
! ---------------------------------------------------------

CONTAINS

!*
! NAME
!     MultiplySfilter
!
! DESCRIPTION
!     Multiply Sfilter types
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION MultiplySfilter(s1,s2) RESULT(res)

USE filter_types
USE filter_functions
USE polynomial_operators

IMPLICIT NONE

! argument types
type(Sfilter), intent(in) :: s1,s2
! Result type
type(Sfilter) :: res

! function definition

res%a%order=s1%a%order+s2%a%order
res%b%order=s1%b%order+s2%b%order

if (s1%wnorm.ne.s2%wnorm) then
  write(*,*)'wnorm discrepancy in MultiplySfilter'
  stop
end if

res=allocate_Sfilter(res%a%order,res%b%order)
res%wnorm=s1%wnorm

! multiply polynomials top and bottom to give the result
res%a=s1%a*s2%a
res%b=s1%b*s2%b

END FUNCTION MultiplySfilter

! ---------------------------------------------------------
!+
! NAME
!     MultiplyZfilter
!
! DESCRIPTION
!     Multiply Zfilter types
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION MultiplyZfilter(z1,z2) RESULT(res)

USE filter_types
USE filter_functions
USE polynomial_operators

IMPLICIT NONE

! argument types
type(Zfilter), intent(in) :: z1,z2
! Result type
type(Zfilter) :: res

! function definition

res%a%order=z1%a%order+z2%a%order
res%b%order=z1%b%order+z2%b%order

if (z1%wnorm.ne.z2%wnorm) then
  write(*,*)'wnorm discrepancy in MultiplyZfilter'
  stop
end if

if (z1%T.ne.z2%T) then
  write(*,*)'Timestep discrepancy in MultiplyZfilter'
  stop
end if

res=allocate_Zfilter(res%a%order,res%b%order)
res%T=z1%T
res%wnorm=z1%wnorm

! multiply polynomials top and bottom to give the result
res%a=z1%a*z2%a
res%b=z1%b*z2%b

END FUNCTION MultiplyZfilter

! ---------------------------------------------------------
!*
! NAME
!     DivideSfilter
!
! DESCRIPTION
!     Divide Sfilter types
!
! HISTORY
!
!     started 26/01/09 CJS
!

FUNCTION DivideSfilter(s1,s2) RESULT(res)

USE filter_types
USE filter_functions
USE polynomial_operators

IMPLICIT NONE

! argument types
type(Sfilter), intent(in) :: s1,s2
! Result type
type(Sfilter) :: res

! function definition

res%a%order=s1%a%order+s2%b%order
res%b%order=s1%b%order+s2%a%order

if (s1%wnorm.ne.s2%wnorm) then
  write(*,*)'wnorm discrepancy in DivideSfilter'
  stop
end if

res=allocate_Sfilter(res%a%order,res%b%order)
res%wnorm=s1%wnorm

! Divide polynomials top and bottom to give the result
res%a=s1%a*s2%b
res%b=s1%b*s2%a

END FUNCTION DivideSfilter

! ---------------------------------------------------------
!+
! NAME
!     DivideZfilter
!
! DESCRIPTION
!     Divide Zfilter types
!
! HISTORY
!
!     started 26/01/09 CJS
!

FUNCTION DivideZfilter(z1,z2) RESULT(res)

USE filter_types
USE filter_functions
USE polynomial_operators

IMPLICIT NONE

! argument types
type(Zfilter), intent(in) :: z1,z2
! Result type
type(Zfilter) :: res

! function definition

res%a%order=z1%a%order+z2%b%order
res%b%order=z1%b%order+z2%a%order

if (z1%wnorm.ne.z2%wnorm) then
  write(*,*)'wnorm discrepancy in DivideZfilter'
  stop
end if

if (z1%T.ne.z2%T) then
  write(*,*)'Timestep discrepancy in DivideZfilter'
  stop
end if

res=allocate_Zfilter(res%a%order,res%b%order)
res%T=z1%T
res%wnorm=z1%wnorm

! Divide polynomials top and bottom to give the result
res%a=z1%a*z2%b
res%b=z1%b*z2%a

END FUNCTION DivideZfilter

! ---------------------------------------------------------

SUBROUTINE AssignSfilter(res,s)
USE filter_types
IMPLICIT NONE

! Result type
type(Sfilter), intent(out) :: res
! argument types
real*8	    , intent(in)  :: s

! SUBROUTINE definition
  res%wnorm=1d0
  res%a%order=0
  allocate (res%a%coeff(0:res%a%order))
  res%b%order=0
  allocate (res%b%coeff(0:res%b%order))
  res%a%coeff(0:res%a%order)=s
  res%b%coeff(0:res%b%order)=1d0

END SUBROUTINE AssignSfilter
! ---------------------------------------------------------

SUBROUTINE AssignZfilter(res,s)
USE filter_types
IMPLICIT NONE

! Result type
type(Zfilter), intent(out) :: res
! argument types
real*8	    , intent(in)  :: s

! SUBROUTINE definition
  res%wnorm=1d0
  res%T=1D0
  res%a%order=0
  allocate (res%a%coeff(0:res%a%order))
  res%b%order=0
  allocate (res%b%coeff(0:res%b%order))
  res%a%coeff(0:res%a%order)=s
  res%b%coeff(0:res%b%order)=1d0

END SUBROUTINE AssignZfilter

! ---------------------------------------------------------

SUBROUTINE AssignSfilter2(res,s)
USE filter_types
IMPLICIT NONE

! Result type
type(Sfilter) , intent(out) :: res
! argument types
type(Sfilter) , intent(in)  :: s

! SUBROUTINE definition

! OLD VERSION, MAY CAUSE AN ERROR IF RES IS ALREADY ALLOCATED
!  if (.NOT. ALLOCATED (res%a%coeff(0:res%a%order))) then
!! this filter not yet allocated so allocate here
!    ALLOCATE (res%a%coeff(0:s%a%order))
!  else if (res%a%order.ne.s%a%order) then
!! this filter has already been allocated but to the wrong order so reallocate  
!    DEALLOCATE (res%a%coeff)
!    ALLOCATE (res%a%coeff(0:s%a%order))
!  end if
!  
!  if (.NOT. ALLOCATED (res%b%coeff(0:res%b%order))) then
!! this filter not yet allocated so allocate here
!    ALLOCATE (res%b%coeff(0:s%b%order))
!  else if (res%b%order.ne.s%b%order) then
!! this filter has already been allocated but to the wrong order so reallocate  
!    DEALLOCATE (res%b%coeff)
!    ALLOCATE (res%b%coeff(0:s%b%order))
!  end if
  
  if (ALLOCATED (res%a%coeff)) DEALLOCATE (res%a%coeff)
  ALLOCATE (res%a%coeff(0:s%a%order))

  if (ALLOCATED (res%b%coeff)) DEALLOCATE (res%b%coeff)
  ALLOCATE (res%b%coeff(0:s%b%order))
  
  
! set values  
  res%wnorm=s%wnorm
  res%a%order=s%a%order
  res%a%coeff(0:res%a%order)=s%a%coeff(0:s%a%order)
  res%b%order=s%b%order
  res%b%coeff(0:res%b%order)=s%b%coeff(0:s%b%order)

END SUBROUTINE AssignSfilter2
! ---------------------------------------------------------

SUBROUTINE AssignZfilter2(res,s)
USE filter_types
IMPLICIT NONE

! Result type
type(Zfilter) , intent(out) :: res
! argument types
type(Zfilter) , intent(in)  :: s

! SUBROUTINE definition
  
!  if (.NOT. ALLOCATED (res%a%coeff(0:res%a%order))) then
  if (.NOT. ALLOCATED (res%a%coeff)) then
! this filter not yet allocated so allocate here
    ALLOCATE (res%a%coeff(0:s%a%order))
  else if (res%a%order.ne.s%a%order) then
! this filter has already been allocated but to the wrong order so reallocate  
    DEALLOCATE (res%a%coeff)
    ALLOCATE (res%a%coeff(0:s%a%order))
  end if
  
!  if (.NOT. ALLOCATED (res%b%coeff(0:res%b%order))) then
  if (.NOT. ALLOCATED (res%b%coeff)) then
! this filter not yet allocated so allocate here
    ALLOCATE (res%b%coeff(0:s%b%order))
  else if (res%b%order.ne.s%b%order) then
! this filter has already been allocated but to the wrong order so reallocate  
    DEALLOCATE (res%b%coeff)
    ALLOCATE (res%b%coeff(0:s%b%order))
  end if
  
! set values  
  res%wnorm=s%wnorm
  res%T=s%T
  res%a%order=s%a%order
  res%a%coeff(0:res%a%order)=s%a%coeff(0:s%a%order)
  res%b%order=s%b%order
  res%b%coeff(0:res%b%order)=s%b%coeff(0:s%b%order)

END SUBROUTINE AssignZfilter2
! ---------------------------------------------------------

SUBROUTINE AssignS_PR_filter(res,s)
USE filter_types
IMPLICIT NONE

! Result type
type(Sfilter_PR) , intent(out) :: res
! argument types
type(Sfilter_PR) , intent(in)  :: s
  
  if (ALLOCATED (res%complex_pole)) DEALLOCATE (res%complex_pole)
  ALLOCATE (res%complex_pole(1:s%order))
  
  if (ALLOCATED (res%poles)) DEALLOCATE (res%poles)
  ALLOCATE (res%poles(1:s%order))
  
  if (ALLOCATED (res%residues)) DEALLOCATE (res%residues)
  ALLOCATE (res%residues(1:s%order))
 
! set values  
  res%wnorm=s%wnorm
  res%order=s%order
  res%C=s%C
  
  res%n_complex_poles=s%n_complex_poles
  res%n_complex_pole_pairs=s%n_complex_pole_pairs
  res%n_real_poles=s%n_real_poles
  res%complex_pole(1:res%order)=s%complex_pole(1:s%order)
  res%poles(1:res%order)=s%poles(1:s%order)

  res%residues(1:res%order)=s%residues(1:s%order) 

END SUBROUTINE AssignS_PR_filter
! ---------------------------------------------------------

SUBROUTINE AssignS_PZ_filter(res,s)
USE filter_types
IMPLICIT NONE

! Result type
type(Sfilter_PZ) , intent(out) :: res
! argument types
type(Sfilter_PZ) , intent(in)  :: s
  
  if (ALLOCATED (res%complex_pole)) DEALLOCATE (res%complex_pole)
  ALLOCATE (res%complex_pole(1:s%order))
  
  if (ALLOCATED (res%poles)) DEALLOCATE (res%poles)
  ALLOCATE (res%poles(1:s%order))
  
  if (ALLOCATED (res%complex_zero)) DEALLOCATE (res%complex_zero)
  ALLOCATE (res%complex_zero(1:s%order))
  
  if (ALLOCATED (res%zeros)) DEALLOCATE (res%zeros)
  ALLOCATE (res%zeros(1:s%order))
  
! set values  
  res%wnorm=s%wnorm
  res%order=s%order
  res%G=s%G
  
  res%n_complex_poles=s%n_complex_poles
  res%n_complex_pole_pairs=s%n_complex_pole_pairs
  res%n_real_poles=s%n_real_poles
  res%complex_pole(1:res%order)=s%complex_pole(1:s%order)
  res%poles(1:res%order)=s%poles(1:s%order)
  
  res%n_complex_zeros=s%n_complex_zeros
  res%n_complex_zero_pairs=s%n_complex_zero_pairs
  res%n_real_zeros=s%n_real_zeros
  res%complex_zero(1:res%order)=s%complex_zero(1:s%order)
  res%zeros(1:res%order)=s%zeros(1:s%order)

END SUBROUTINE AssignS_PZ_filter

END MODULE filter_operators



