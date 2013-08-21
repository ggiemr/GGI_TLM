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
!     Vars
!
! Description
!     type definitions and global variables for filter representations
!
! Comments:
!      
!
! History
!
!     started 22/01/09 CJS
!

MODULE filter_types

USE polynomial_types

TYPE::Sfilter

  REAL*8  :: wnorm
  TYPE(polynomial) :: a
  TYPE(polynomial) :: b
  
END TYPE Sfilter

TYPE::Zfilter

! Note the polynomials are in terms of z^(-1) 

  REAL*8  :: wnorm
  REAL*8  :: T               ! Note this is a normalised timestep
  TYPE(polynomial) :: a
  TYPE(polynomial) :: b
  
END TYPE Zfilter

TYPE::Zfilter_response  ! response of filter f(z)=(a(z)/b(z)) g(z)

  INTEGER	     :: order
  REAL*8	     :: f
  REAL*8,allocatable :: w(:)
  
END TYPE Zfilter_response

TYPE::Sfilter_PZ

  REAL*8  			:: wnorm
  INTEGER			:: order
  REAL*8			:: G
  INTEGER			:: n_complex_poles
  INTEGER			:: n_complex_pole_pairs
  INTEGER			:: n_real_poles
  LOGICAL,allocatable		:: complex_pole(:)
  COMPLEX*16,allocatable	:: poles(:)
  INTEGER			:: n_complex_zeros
  INTEGER			:: n_complex_zero_pairs
  INTEGER			:: n_real_zeros
  LOGICAL,allocatable		:: complex_zero(:)
  COMPLEX*16,allocatable	:: zeros(:)
  
END TYPE Sfilter_PZ

TYPE::Sfilter_PR

  REAL*8  			:: wnorm
  INTEGER			:: order
  REAL*8			:: C
  INTEGER			:: n_complex_poles
  INTEGER			:: n_complex_pole_pairs
  INTEGER			:: n_real_poles
  LOGICAL,allocatable		:: complex_pole(:)
  COMPLEX*16,allocatable	:: poles(:)
  COMPLEX*16,allocatable	:: residues(:)
  
END TYPE Sfilter_PR

END MODULE filter_types

