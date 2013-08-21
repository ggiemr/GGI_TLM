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
MODULE pul_wire_data

IMPLICIT NONE

TYPE::pul_wire_geometry

  integer :: nterms
  integer :: npoints

  real*8  :: xc
  real*8  :: yc
  real*8  :: rw
  real*8  :: ri
  real*8  :: epsr
  
  real*8,allocatable  :: xp(:)
  real*8,allocatable  :: yp(:)
  real*8,allocatable  :: rp(:)
  real*8,allocatable  :: tp(:)
  
  real*8		:: alpha
  real*8,allocatable  	:: a(:)
  real*8,allocatable  	:: b(:)
  integer 		:: matrix_position
  integer 		:: matrix_position2
  
  real*8,allocatable  :: xdp(:)
  real*8,allocatable  :: ydp(:)
  real*8,allocatable  :: rdp(:)
  real*8,allocatable  :: tdp(:)
  
  real*8,allocatable  :: nxdp(:)
  real*8,allocatable  :: nydp(:)
  
  real*8		:: alpha2
  real*8,allocatable  	:: a2(:)
  real*8,allocatable  	:: b2(:)
  
END TYPE pul_wire_geometry


END MODULE pul_wire_data


MODULE pul_data

USE pul_wire_data

IMPLICIT NONE

  integer :: pul_nwires_in
  integer :: pul_nwires
  
  real*8  :: pul_dl
  
  real*8  :: pul_xmin,pul_xmax
  integer :: pul_nxp
  
  real*8  :: pul_ymin,pul_ymax
  integer :: pul_nyp
  
  integer :: pul_tot_nterms
  integer :: pul_tot_nterms_dielectric
  integer :: pul_tot_npoints
  integer :: pul_tot_npoints_dielectric
  integer :: pul_return_conductor
  
  integer :: pul_LC_matdim
  integer :: pul_D_matdim
  
  logical :: pul_include_dielectric
  logical :: pul_include_TLM_return
  logical :: pul_op_flag
  
  TYPE(pul_wire_geometry),allocatable  :: pul_wire_spec(:)

  real*8,allocatable  :: pul_D(:,:)
  real*8,allocatable  :: pul_B(:,:)
  
  real*8,allocatable  :: pul_phi(:)
  real*8,allocatable  :: pul_alpha(:)
  
  real*8,allocatable  :: pul_L(:,:)
  real*8,allocatable  :: pul_C(:,:)
  real*8,allocatable  :: pul_Cg(:,:)

END MODULE pul_data
