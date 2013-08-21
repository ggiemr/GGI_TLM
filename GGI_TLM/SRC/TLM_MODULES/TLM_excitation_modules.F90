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
! MODULE TLM_excitation
!
! NAME
!     MODULE TLM_excitation
!
! DESCRIPTION
!     excitation data relating to the TLM solution
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
MODULE TLM_excitation

USE geometry_types
!USE mode_field ! no longer used

IMPLICIT NONE

TYPE::excitation_function_type

  integer	:: type
  integer	:: n_parameters
  real*8	:: parameters(6) 
  real*8,allocatable	:: value(:)
  real*8,allocatable	:: value_face(:)
 
END TYPE excitation_function_type

TYPE::excitation_point_type

  integer		:: excitation_function_number
  integer		:: point_number
  type(xyz)		:: point
  type(cell_point)	:: cell_point
  integer		:: rank
  integer		:: cell_excitation_field_number
  integer		:: face_excitation_field_number
  integer		:: field_component
  integer		:: source_type
   
END TYPE excitation_point_type

TYPE::excitation_surface_type

  integer			:: excitation_function_number
  integer			:: surface_number
  logical 			:: excitation_on_outward_normal
  integer 			:: number_of_faces
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: face_excitation_field_number_list(:)
  integer			:: field_component
  integer			:: source_type
   
END TYPE excitation_surface_type

TYPE::huygens_surface_type

  integer		:: surface_number
  logical 		:: excitation_on_outward_normal
  integer		:: excitation_function_number
  real*8		:: Ktheta
  real*8		:: Kphi
  real*8		:: Etheta
  real*8		:: Ephi

  logical       	:: outer_surface_flag
  integer       	:: n_surface_patches
  real*8        	:: Ei(3),Hi(3)
  real*8        	:: Ki(3)
  real*8        	:: xmin,xmax,ymin,ymax,zmin,zmax
  real*8        	:: offset_min
  real*8,allocatable	:: offset(:)
  integer,allocatable	:: cx(:)
  integer,allocatable	:: cy(:)
  integer,allocatable	:: cz(:)
  integer,allocatable	:: face(:)
  integer,allocatable	:: nx(:)
  integer,allocatable	:: ny(:)
  integer,allocatable	:: nz(:)  
  integer,allocatable	:: face_excitation_field_number_list(:)

END TYPE huygens_surface_type

TYPE::excitation_mode_type

  INTEGER	    	:: excitation_function_number
  INTEGER	    	:: surface_number
  logical 		:: excitation_on_outward_normal
  integer		:: field_component
  integer		:: source_type
  CHARACTER*256		:: mode_file_name 
  
  INTEGER			::  n_mode_samples
  
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: face_excitation_field_number_list(:)
  real*8,allocatable		:: mode_field(:)

  INTEGER		:: xcol
  INTEGER		:: ycol
  INTEGER		:: zcol
  INTEGER		:: mode_col
  
END TYPE excitation_mode_type
  
  integer,parameter	:: excitation_function_type_impulse=1
  integer,parameter	:: excitation_function_type_gaussian=2
  integer,parameter	:: excitation_function_type_gaussian_step=3
  integer,parameter	:: excitation_function_type_step=4
  integer,parameter	:: excitation_function_type_sinusoid=5
  integer,parameter	:: excitation_function_type_gaussian_sinusoid=6
  integer,parameter	:: excitation_function_type_gaussian_step_sinusoid=7
  integer,parameter	:: excitation_function_type_double_exponential=8
  
  integer,parameter	:: source_type_hard=1
  
  integer,parameter	:: source_type_soft=2

  integer		:: total_number_excitation_cells
  real*8,allocatable	:: cell_excitation_field(:,:)

  integer		:: total_number_excitation_faces
  real*8,allocatable	:: face_excitation_field(:,:,:)
  
  integer				     	:: n_excitation_functions
  type(excitation_function_type),allocatable 	:: excitation_functions(:)
  
  integer				     	:: n_excitation_points
  type(excitation_point_type),allocatable 	:: excitation_points(:)
  
  integer				     	:: n_excitation_surfaces
  type(excitation_surface_type),allocatable 	:: excitation_surfaces(:)
  
  integer					:: n_huygens_surfaces
  type(huygens_surface_type) 			:: huygens_surface
  
  integer					:: n_excitation_modes
  type(excitation_mode_type),allocatable 	:: excitation_mode_list(:)


END MODULE TLM_excitation
