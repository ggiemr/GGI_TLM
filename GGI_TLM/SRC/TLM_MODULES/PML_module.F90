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
! NAME
!     MODULE PML_module
!
! DESCRIPTION
!     PML information
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/10/2019 CJS
!
!
MODULE PML_module

USE geometry_types

IMPLICIT NONE

integer :: n_pml_volumes

logical :: pml_xmin_flag=.FALSE.
logical :: pml_xmax_flag=.FALSE.
logical :: pml_ymin_flag=.FALSE.
logical :: pml_ymax_flag=.FALSE.
logical :: pml_zmin_flag=.FALSE.
logical :: pml_zmax_flag=.FALSE.

integer :: pml_volume_to_face(1:6)

integer,parameter :: pml_face_xmin=1
integer,parameter :: pml_face_xmax=2
integer,parameter :: pml_face_ymin=3
integer,parameter :: pml_face_ymax=4
integer,parameter :: pml_face_zmin=5
integer,parameter :: pml_face_zmax=6

real*8  :: pml_txmin,pml_txmax,pml_tymin,pml_tymax,pml_tzmin,pml_tzmax
  
type(volume_type),allocatable 	:: pml_volumes(:)

real*8  :: pml_r
   
END MODULE PML_module
