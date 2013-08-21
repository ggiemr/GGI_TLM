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
! NAME
!     MODULE mode_stir
!
! DESCRIPTION
!     data relating to mode stir boundary conditions
!	
!     
! COMMENTS
!  
!
!
! HISTORY
!
!     started 10/01/13 CJS
!
!
MODULE mode_stir

USE geometry_types

IMPLICIT NONE

TYPE::mode_stir_surface_type
  
  REAL*8   		:: R_ms
  REAL*8   		:: impulse_amplitude
  integer 		:: number_of_surfaces
  integer,allocatable	:: surface_list(:)
  integer,allocatable	:: surface_orientation_list(:)
  
  integer 		:: number_of_faces
  integer 		:: number_of_ports
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: port1(:)
  integer,allocatable		:: port2(:)
  real*8,allocatable		:: port_voltage_list(:)
  
  integer,allocatable			:: number_of_ports_rank(:)
  
  integer 				:: total_number_of_mode_stir_faces
  integer 				:: total_number_of_mode_stir_ports
  integer,allocatable			:: mode_stir_voltage_pairing(:)
  integer,allocatable			:: sign(:)
  
END TYPE mode_stir_surface_type
  
integer					 :: n_mode_stir_surfaces
type(mode_stir_surface_type),allocatable :: mode_stir_surface_list(:)

END MODULE mode_stir
