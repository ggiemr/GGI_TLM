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
! MODULE post_process
!
! NAME
!     MODULE post_process
!
! DESCRIPTION
!     post processing data structures and data
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/08/2012 CJS
!
!
MODULE post_process

IMPLICIT NONE
  
TYPE::time_domain_data

  integer		:: n_timesteps
  real*8,allocatable	:: time(:)
  real*8,allocatable	:: value(:)
 
END TYPE time_domain_data

  
TYPE::frequency_domain_data

  integer			:: n_frequencies
  real*8,allocatable		:: frequency(:)
  complex*16,allocatable	:: value(:)
  real*8,allocatable		:: magnitude(:)
  real*8,allocatable		:: phase(:)
  real*8,allocatable		:: dB(:)
 
END TYPE frequency_domain_data
  
TYPE::surface_animation_data

  integer			:: n_points
  integer			:: n_quads
  integer			:: n_frames
  
  real*8,allocatable 		:: points(:,:)

  integer,allocatable 		:: quads(:,:)
  
  real*8,allocatable 		:: frame_data(:,:) 
   
  complex*16,allocatable 	:: complex_data(:)  
  complex*16,allocatable 	:: magnitude_data(:)  
  integer,allocatable 		:: cx(:)
  integer,allocatable 		:: cy(:)
  integer,allocatable 		:: cz(:)
  
  real*8 			:: max_data
  real*8 			:: min_data

END TYPE surface_animation_data

integer					:: n_functions_of_time
TYPE(time_domain_data),allocatable	:: function_of_time(:)

integer					:: n_functions_of_frequency
TYPE(frequency_domain_data),allocatable	:: function_of_frequency(:)
  
END MODULE post_process

