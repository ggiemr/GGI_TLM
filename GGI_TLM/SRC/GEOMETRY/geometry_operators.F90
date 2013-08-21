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
! Name
!     Operators -  overloading operators for geometry types
!
! Description
!     This module defines arithmetic operations for geometry types
!
! Public
!     The following are defined in this module
!
!     +
!     -
!     =
!
! History
!
!     started 11/09/12 CJS
!
MODULE geometry_operators

IMPLICIT NONE

PRIVATE  ! everything private unless stated explicitly

PUBLIC :: OPERATOR(+), OPERATOR(-), ASSIGNMENT(=)

! OVERLOAD + to operate on type(xyz) 

INTERFACE OPERATOR(+)
  MODULE PROCEDURE Add_xyz
END INTERFACE

! OVERLOAD - to operate on type(xyz)

INTERFACE OPERATOR(-)
  MODULE PROCEDURE Subtract_xyz
END INTERFACE

! OVERLOAD = to assign type(xyz)

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE Assign_xyz,Assign_xyz_line,Assign_xyz_triangle,Assign_xyz_tet,Assign_ijk,Assign_cell_point
END INTERFACE


! ---------------------------------------------------------

CONTAINS

!+
! NAME
!     Add_xyz
!
! DESCRIPTION
!     Add xyz types
!
! HISTORY
!
!     started 11/09/12 CJS
!

FUNCTION Add_xyz(xyz1,xyz2) RESULT(res)

USE geometry_types

IMPLICIT NONE

! argument types
type(xyz), intent(in) :: xyz1,xyz2
! Result type
type(xyz) :: res

! START

res%x=xyz1%x+xyz2%x
res%y=xyz1%y+xyz2%y
res%z=xyz1%z+xyz2%z

END FUNCTION add_xyz

! ---------------------------------------------------------

!+
! NAME
!     subtract_xyz
!
! DESCRIPTION
!     subtract xyz types
!
! HISTORY
!
!     started 11/09/12 CJS
!

FUNCTION subtract_xyz(xyz1,xyz2) RESULT(res)

USE geometry_types

IMPLICIT NONE

! argument types
type(xyz), intent(in) :: xyz1,xyz2
! Result type
type(xyz) :: res

! START

res%x=xyz1%x-xyz2%x
res%y=xyz1%y-xyz2%y
res%z=xyz1%z-xyz2%z

END FUNCTION subtract_xyz

! 
! ---------------------------------------------------------
! NAME
!     assign_xyz
!
! DESCRIPTION
!     assign xyz types
!
! COMMENTS
!     
! HISTORY
!
!     started 11/09/12 CJS
!

SUBROUTINE Assign_xyz(res,b)

USE geometry_types

IMPLICIT NONE

! Result type
type(xyz), intent(out):: res
! argument types
type(xyz), intent(in) :: b

! START

  res%x=b%x
  res%y=b%y
  res%z=b%z

END SUBROUTINE Assign_xyz

! 
! ---------------------------------------------------------
! NAME
!     assign_xyz_line
!
! DESCRIPTION
!     assign xyz_line types
!
! COMMENTS
!     
! HISTORY
!
!     started 13/09/12 CJS
!

SUBROUTINE Assign_xyz_line(res,b)

USE geometry_types

IMPLICIT NONE

! Result type
type(xyz_line), intent(out):: res
! argument types
type(xyz_line), intent(in) :: b

! START

  res%end(1)%x=b%end(1)%x
  res%end(1)%y=b%end(1)%y
  res%end(1)%z=b%end(1)%z
  res%end(2)%x=b%end(2)%x
  res%end(2)%y=b%end(2)%y
  res%end(2)%z=b%end(2)%z

END SUBROUTINE Assign_xyz_line

! 
! ---------------------------------------------------------
! NAME
!     assign_xyz_triangle
!
! DESCRIPTION
!     assign xyz_triangle types
!
! COMMENTS
!     
! HISTORY
!
!     started 13/09/12 CJS
!

SUBROUTINE Assign_xyz_triangle(res,b)

USE geometry_types

IMPLICIT NONE

! Result type
type(xyz_triangle), intent(out):: res
! argument types
type(xyz_triangle), intent(in) :: b

! START

  res%vertex(1)%x=b%vertex(1)%x
  res%vertex(1)%y=b%vertex(1)%y
  res%vertex(1)%z=b%vertex(1)%z
  res%vertex(2)%x=b%vertex(2)%x
  res%vertex(2)%y=b%vertex(2)%y
  res%vertex(2)%z=b%vertex(2)%z
  res%vertex(3)%x=b%vertex(3)%x
  res%vertex(3)%y=b%vertex(3)%y
  res%vertex(3)%z=b%vertex(3)%z

END SUBROUTINE Assign_xyz_triangle

! 
! ---------------------------------------------------------
! NAME
!     assign_xyz_tet
!
! DESCRIPTION
!     assign xyz_tet types
!
! COMMENTS
!     
! HISTORY
!
!     started 13/09/12 CJS
!

SUBROUTINE Assign_xyz_tet(res,b)

USE geometry_types

IMPLICIT NONE

! Result type
type(xyz_tet), intent(out):: res
! argument types
type(xyz_tet), intent(in) :: b

! START

  res%vertex(1)%x=b%vertex(1)%x
  res%vertex(1)%y=b%vertex(1)%y
  res%vertex(1)%z=b%vertex(1)%z
  res%vertex(2)%x=b%vertex(2)%x
  res%vertex(2)%y=b%vertex(2)%y
  res%vertex(2)%z=b%vertex(2)%z
  res%vertex(3)%x=b%vertex(3)%x
  res%vertex(3)%y=b%vertex(3)%y
  res%vertex(3)%z=b%vertex(3)%z
  res%vertex(4)%x=b%vertex(4)%x
  res%vertex(4)%y=b%vertex(4)%y
  res%vertex(4)%z=b%vertex(4)%z

END SUBROUTINE Assign_xyz_tet

! 
! ---------------------------------------------------------
! NAME
!     assign_ijk
!
! DESCRIPTION
!     assign ijk types
!
! COMMENTS
!     
! HISTORY
!
!     started 13/09/12 CJS
!

SUBROUTINE Assign_ijk(res,b)

USE geometry_types

IMPLICIT NONE

! Result type
type(ijk), intent(out):: res
! argument types
type(ijk), intent(in) :: b

! START

  res%i=b%i
  res%j=b%j
  res%k=b%k

END SUBROUTINE Assign_ijk
! 
! ---------------------------------------------------------
! NAME
!     assign_cell_point
!
! DESCRIPTION
!     assign cell_point types
!
! COMMENTS
!     
! HISTORY
!
!     started 13/09/12 CJS
!

SUBROUTINE Assign_cell_point(res,b)

USE geometry_types

IMPLICIT NONE

! Result type
type(cell_point), intent(out):: res
! argument types
type(cell_point), intent(in) :: b

! START

  res%cell%i=b%cell%i
  res%cell%j=b%cell%j
  res%cell%k=b%cell%k
  res%point =b%point

END SUBROUTINE Assign_cell_point

END MODULE  geometry_operators

