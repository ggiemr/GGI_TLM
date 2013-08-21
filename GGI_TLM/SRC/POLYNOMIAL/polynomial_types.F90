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
! Name
!    polynomial_types
!
! Description
!     type definitions for polynomials
!
! Comments:
!      
!
! History
!
!     started 23/01/09 CJS
!

MODULE polynomial_types

TYPE::polynomial
  INTEGER	     :: order
  REAL*8,allocatable :: coeff(:)
END TYPE polynomial

TYPE::complex_polynomial
  INTEGER	     :: order
  COMPLEX*16,allocatable :: coeff(:)
END TYPE complex_polynomial

END MODULE polynomial_types

