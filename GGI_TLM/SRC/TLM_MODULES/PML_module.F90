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

! Notation change for PML to make things consistent with the paper
! PML port numbering               ! GGI_TLM port numbering

integer,parameter :: Vxny=1   !  Vy_xmin=1
integer,parameter :: Vxnz=2   !  Vz_xmin=2
integer,parameter :: Vxpy=3   !  Vy_xmax=3
integer,parameter :: Vxpz=4   !  Vz_xmax=4
integer,parameter :: Vynx=5   !  Vx_ymin=5
integer,parameter :: Vynz=6   !  Vz_ymin=6
integer,parameter :: Vypx=7   !  Vx_ymax=7
integer,parameter :: Vypz=8   !  Vz_ymax=8
integer,parameter :: Vznx=9   !  Vx_zmin=9
integer,parameter :: Vzny=10  !  Vy_zmin=10
integer,parameter :: Vzpx=11  !  Vx_zmax=11
integer,parameter :: Vzpy=12  !  Vy_zmax=12

real*8  :: pml_txmin,pml_txmax,pml_tymin,pml_tymax,pml_tzmin,pml_tzmax
  
type(volume_type),allocatable 	:: pml_volumes(:)

! PML estimated reflection coefficient
real*8  :: pml_r_xmin,pml_r_xmax
real*8  :: pml_r_ymin,pml_r_ymax
real*8  :: pml_r_zmin,pml_r_zmax

! PML order
integer :: pml_order

! conductivity of first PML layer in x, y and z directions
real*8  :: pml_s0_xmin,pml_s0_xmax
real*8  :: pml_s0_ymin,pml_s0_ymax
real*8  :: pml_s0_zmin,pml_s0_zmax

logical ::   PML_material_intersection_flag
logical ::   PML_cable_intersection_flag
logical ::   PML_excitation_intersection_flag

integer :: pml_xmin,pml_xmax,npml_x
integer :: pml_ymin,pml_ymax,npml_y
integer :: pml_zmin,pml_zmax,npml_z

! array to transform from d_x,d_y,d_z cell to position in the PML parameters array

integer,allocatable :: PML_dxdydz_to_parameter_array(:,:,:)

! Parameters for each type of PML cell i.e. for each sx, sy, sz combination
TYPE::PML_material_type

  integer :: d_x,d_y,d_z  ! depth into PML in x, y and z for this cell
  real*8  :: ax,ay,az     ! alpha constants of each cell type
  real*8  :: sx,sy,sz     ! conductivity of each cell type

END TYPE PML_material_type

! data for each PML cell i.e. voltage and current values to be saved for each cell
TYPE::PML_data_type

  integer :: PML_parameter_array_pos ! position in PML_parameters array to find the correct constants for this cell
  
  real*8  :: Vxt,Vyt,Vzt
  real*8  :: Ix,Iy,Iz
  real*8  :: Vi(12)
  real*8  :: Vyxt,Vzxt,Vxyt,Vzyt,Vyzt,Vxzt

END TYPE PML_data_type

integer :: PML_array_size
integer :: PML_n_parameters

type(PML_material_type),allocatable :: PML_parameters(:)

integer	:: total_number_of_PML_cells

type(PML_data_type),allocatable :: PML_cell_data(:)

END MODULE PML_module
