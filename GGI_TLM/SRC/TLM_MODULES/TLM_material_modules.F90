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
! MODULE TLM_volume_materials
! MODULE TLM_surface_materials
!
! NAME
!     MODULE TLM_volume_materials
!
! DESCRIPTION
!     volume_materials data relating to the TLM solution
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/09/2012 CJS
!
!
MODULE TLM_volume_materials


USE filter_types

IMPLICIT NONE

TYPE::volume_material_type

  character(len=256)	:: name
  
  integer		:: type
  integer		:: n_volumes
  integer,allocatable	:: volume_list(:)
  
  REAL*8 		:: fmin,fmax
  
  logical		:: eps_filter_exists
  TYPE(Sfilter)		:: eps_S
  TYPE(Sfilter)		:: Zs_eps_S
  TYPE(Zfilter)		:: Zs_eps_Z
  REAL*8	       	:: Zs_eps_f
  REAL*8	       	:: Ys_eps_f
  
  logical		:: mu_filter_exists
  TYPE(Sfilter)		:: mu_S  
  TYPE(Sfilter)		:: Zs_mu_S  
  TYPE(Zfilter)		:: Zs_mu_Z  
  REAL*8	       	:: Zs_mu_f
  
  REAL*8		:: sigma_e
  REAL*8		:: sigma_m
  REAL*8	       	:: Ge
  REAL*8	       	:: Rm

END TYPE volume_material_type

  integer,parameter	:: volume_material_type_PEC=1

  integer,parameter	:: volume_material_type_PMC=2

  integer,parameter	:: volume_material_type_DISPERSIVE=3

  integer				  :: n_volume_materials
  type(volume_material_type),allocatable :: volume_material_list(:)
  
  integer 	:: n_volume_material_cells
  integer 	:: volume_material_storage
  
  type(Zfilter_response),allocatable	:: volume_material_Zs_eps_filter_data(:) 
  type(Zfilter_response),allocatable	:: volume_material_Zs_mu_filter_data(:) 
  
END MODULE TLM_volume_materials
!
! NAME
!     MODULE TLM_surface_materials
!
! DESCRIPTION
!     surface_materials data relating to the TLM solution
!     
! COMMENTS
!     
!
! HISTORY
!
!    started 10/08/2012 CJS
!    2/12/2013 		CJS: Implement anisotropic impedance boundary conditions
!
!
MODULE TLM_surface_materials

USE filter_types

IMPLICIT NONE

TYPE::surface_material_type

  character(len=256)	:: name

  integer		:: type
  integer		:: n_surfaces
  integer,allocatable	:: surface_list(:)
  integer,allocatable	:: surface_orientation_list(:)
  
  REAL*8 		:: fmin,fmax
  
  REAL*8 		:: Z11_f(3),Z12_f(3),Z21_f(3),Z22_f(3)
  TYPE(Sfilter)		:: Z11_S(3),Z12_S(3),Z21_S(3),Z22_S(3)
  TYPE(Zfilter)		:: Z11_Z(3),Z12_Z(3),Z21_Z(3),Z22_Z(3)
  
  real*8		:: Diode_Is
  real*8		:: Diode_nVt
  real*8		:: Diode_Rs
  real*8		:: Diode_Cj
  character*2		:: Diode_direction
  integer		:: diode_sign
  
  REAL*8 		:: Diode_Cj_f
  TYPE(Sfilter)		:: Diode_Cj_S
  TYPE(Zfilter)		:: Diode_Cj_Z
   
END TYPE surface_material_type

  integer,parameter	:: surface_material_type_PEC=1

  integer,parameter	:: surface_material_type_PMC=2

  integer,parameter	:: surface_material_type_DISPERSIVE=3

  integer,parameter	:: surface_material_type_FREE_SPACE=4

  integer,parameter	:: surface_material_type_ANISOTROPIC_DISPERSIVE=5
  
  integer,parameter	:: surface_material_type_DIODE=6

  integer				  :: n_surface_materials
  type(surface_material_type),allocatable :: surface_material_list(:)
  
  integer 	:: n_surface_material_faces
  integer 	:: surface_material_storage
  
  type(Zfilter_response),allocatable	:: surface_material_Z11_filter_data(:) 
  type(Zfilter_response),allocatable	:: surface_material_Z12_filter_data(:) 
  type(Zfilter_response),allocatable	:: surface_material_Z21_filter_data(:) 
  type(Zfilter_response),allocatable	:: surface_material_Z22_filter_data(:) 
  
  integer 	:: n_diode_faces
  integer 	:: surface_diode_storage
  type(Zfilter_response),allocatable	:: Diode_Cj_filter_data(:) 

END MODULE TLM_surface_materials
