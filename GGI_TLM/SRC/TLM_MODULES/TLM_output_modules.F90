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
! MODULE TLM_output
!
! NAME
!     MODULE TLM_output
!
! DESCRIPTION
!     output data relating to the TLM solution
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
MODULE TLM_output

USE geometry_types
!USE mode_field ! no longer used

IMPLICIT NONE

TYPE::output_point_type

  integer 			:: point_number
  type(xyz)			:: point
  type(cell_point) 		:: cell_point
  integer			:: rank
  integer			:: field_component
  integer			:: cell_output_field_number
  integer			:: face_output_field_number

  logical			:: specified_timestep_information
  integer			:: first_timestep
  integer			:: last_timestep
  integer			:: timestep_interval
  integer			:: number_of_output_timesteps
  
  logical			:: specified_time_information
  real*8			:: first_time
  real*8			:: last_time
  real*8			:: time_interval
  
  real*8			:: time
  real*8			:: value
   
END TYPE output_point_type

TYPE::output_surface_type

  integer 			:: surface_number
  logical 			:: output_on_outward_normal
  integer 			:: number_of_faces
  type(cell_point),allocatable	:: face_list(:)
  integer			:: field_component
  integer,allocatable		:: face_output_field_number_list(:)
  real*8,allocatable		:: value(:)
  
  logical			:: specified_timestep_information
  integer			:: first_timestep
  integer			:: last_timestep
  integer			:: timestep_interval
  integer			:: number_of_output_timesteps
  
  logical			:: specified_time_information
  real*8			:: first_time
  real*8			:: last_time
  real*8			:: time_interval
  
  integer			:: frame_number
   
END TYPE output_surface_type

TYPE::output_volume_type

  integer 			:: volume_number
  integer 			:: number_of_cells
  type(ijk),allocatable		:: cell_list(:)
  integer,allocatable		:: cell_output_field_number_list(:)
  integer			:: field_component
  real*8,allocatable		:: value(:)
  
  logical			:: specified_timestep_information
  integer			:: first_timestep
  integer			:: last_timestep
  integer			:: timestep_interval
  integer			:: number_of_output_timesteps
  
  logical			:: specified_time_information
  real*8			:: first_time
  real*8			:: last_time
  real*8			:: time_interval
  
  integer			:: frame_number
   
END TYPE output_volume_type

TYPE::output_volume_average_type

  integer 			:: volume_number
  integer 			:: number_of_cells
  type(ijk),allocatable		:: cell_list(:)
  integer,allocatable		:: cell_output_field_number_list(:)
  integer			:: field_component
  real*8			:: value
  
  logical			:: specified_timestep_information
  integer			:: first_timestep
  integer			:: last_timestep
  integer			:: timestep_interval
  integer			:: number_of_output_timesteps
  
  logical			:: specified_time_information
  real*8			:: first_time
  real*8			:: last_time
  real*8			:: time_interval
  
  integer			:: frame_number
   
END TYPE output_volume_average_type

TYPE::far_field_surface_type

  integer 			:: surface_number
  logical 			:: output_on_outward_normal
  real*8  			:: frequency
  real*8  			:: theta_min,theta_max,theta_step
  real*8  			:: phi_min,phi_max,phi_step

  integer 			:: number_of_faces
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: face_output_field_number_list(:)

  complex*16,allocatable	:: J(:,:),M(:,:)
   
END TYPE far_field_surface_type

TYPE::frequency_output_surface_type

  integer 			:: surface_number
  logical 			:: output_on_outward_normal
  real*8  			:: frequency
  
  integer 			:: number_of_faces
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: face_output_field_number_list(:)
  
  integer			:: field_component
  complex*16,allocatable	:: value(:)
   
END TYPE frequency_output_surface_type

TYPE::frequency_output_volume_type

  integer 			:: volume_number
  real*8  			:: frequency
  
  integer 			:: number_of_cells
  type(ijk),allocatable		:: cell_list(:)
  integer,allocatable		:: cell_output_field_number_list(:)
  
  integer			:: field_component
  complex*16,allocatable	:: value(:)
   
END TYPE frequency_output_volume_type

TYPE::frequency_domain_power_surface_type

  integer 			:: surface_number
  logical 			:: output_on_outward_normal
  real*8  			:: fmin
  real*8  			:: fmax
  real*8  			:: fstep
  integer 			:: n_frequencies
  
  integer 			:: number_of_faces
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: face_output_field_number_list(:)
  
  complex*16,allocatable	:: E1(:,:)
  complex*16,allocatable	:: E2(:,:)
  complex*16,allocatable	:: H1(:,:)
  complex*16,allocatable	:: H2(:,:)
  complex*16,allocatable	:: Power(:)
   
END TYPE frequency_domain_power_surface_type

TYPE::RCS_surface_type

  INTEGER	    	:: surface_number
  LOGICAL	    	:: output_on_outward_normal
  REAL*8		:: fmin,fmax,fstep
  REAL*8		:: theta
  REAL*8		:: phi
  
  integer 			:: number_of_faces
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: face_output_field_number_list(:)

  INTEGER   		:: n_far_field_points
  INTEGER   		:: far_field_point_offset
  REAL*8  		:: Kmrcs(3)
  REAL*8  		:: Vtheta(3)
  REAL*8  		:: Vphi(3)
  REAL*8, allocatable   :: Etheta(:)
  REAL*8, allocatable   :: Ephi(:)  
  REAL*8, allocatable   :: Htheta(:)
  REAL*8, allocatable   :: Hphi(:)  
     
END TYPE RCS_surface_type

TYPE::SAR_volume_type

  INTEGER 		:: material_number
  REAL*8		:: frequency
  REAL*8		:: density
  REAL*8		:: conductivity
  REAL*8		:: cell_volume
  REAL*8		:: mass
  INTEGER 		:: number_of_cells
  type(ijk),allocatable	:: cell_list(:)
  integer,allocatable	:: cell_output_field_number_list(:)
  COMPLEX*16,ALLOCATABLE :: Ex(:),Ey(:),Ez(:)
  REAL*8		:: SAR
  
END TYPE SAR_volume_type

TYPE::output_mode_type

  INTEGER	    	:: surface_number
  INTEGER		:: field_component
  CHARACTER*256		:: mode_file_name 
  logical 		:: output_on_outward_normal
  
  INTEGER			::  n_mode_samples
  
  type(cell_point),allocatable	:: face_list(:)
  integer,allocatable		:: face_output_field_number_list(:)
  real*8,allocatable		:: mode_field(:)

  INTEGER		:: xcol
  INTEGER		:: ycol
  INTEGER		:: zcol
  INTEGER		:: mode_col
  
END TYPE output_mode_type
  
integer				    	:: n_output_points
type(output_point_type),allocatable    	:: output_points(:)

integer				    	:: n_output_surfaces
type(output_surface_type),allocatable   :: output_surfaces(:)

integer				    	:: n_output_volumes
type(output_volume_type),allocatable    :: output_volumes(:)

integer				    		:: n_output_volume_averages
type(output_volume_average_type),allocatable    :: output_volume_averages(:)

integer			:: total_number_output_cells
real*8,allocatable	:: cell_output_field(:,:)

integer			:: total_number_output_faces
real*8,allocatable	:: face_output_field(:,:,:)
 
integer					:: n_far_field_surfaces
type(far_field_surface_type)    	:: far_field_surface
 
integer							:: n_frequency_output_surfaces
type(frequency_output_surface_type),allocatable    	:: frequency_output_surface(:)
 
integer							:: n_frequency_output_volumes
type(frequency_output_volume_type),allocatable    	:: frequency_output_volume(:)
 
integer							:: n_frequency_domain_power_surfaces
type(frequency_domain_power_surface_type),allocatable   :: frequency_domain_power_surface(:)

integer					:: n_rcs_surfaces
type(RCS_surface_type)    	:: rcs_surface
integer, parameter 			:: n_rcs_points_per_cell=2

integer				    	:: n_SAR_volumes
type(SAR_volume_type),allocatable    	:: SAR_volume_list(:)
  
integer					:: n_output_modes
type(output_mode_type),allocatable 	:: output_mode_list(:)
 
END MODULE TLM_output
