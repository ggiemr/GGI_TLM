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
! SUBROUTINE reset_packet_data
!
! NAME
!     reset_packet_data
!
! DESCRIPTION
!     reset all packet data
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 12/10/2012 CJS
!
!
SUBROUTINE reset_packet_data

USE TLM_general
USE TLM_volume_materials
USE TLM_surface_materials
USE TLM_excitation
USE TLM_output
USE Cables
USE geometry

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: reset_packet_data',0,output_to_screen_flag)
  
  n_volume_materials=0
  n_surface_materials=0
  n_excitation_functions=0
  n_excitation_points=0
  n_excitation_surfaces=0
  n_huygens_surfaces=0
  n_output_points=0
  n_output_surfaces=0
  n_volumes=0
  n_surfaces=0
  n_lines=0
  n_points=0
  n_cable_geometries=0
  n_cables=0
  n_cable_outputs=0
  n_cable_junctions=0
  n_far_field_surfaces=0
  n_frequency_output_surfaces=0
  n_frequency_output_volumes=0
  n_frequency_domain_power_surfaces=0
  n_output_volumes=0
  n_output_volume_averages=0
  
  bicubic_warp_flag=.FALSE.
  frequency_scale_flag=.FALSE.
  frequency_scale=1d0
  
  wrap_x=.FALSE.
  wrap_y=.FALSE.
  wrap_z=.FALSE.
  
  simulation_time=0d0

! Set default soultion parameters

  new_mesh_generation=.FALSE.
!  new_mesh_generation=.TRUE.

  Cable_LC_Correction_type=LC_correction_type_geometry_scale
  
  number_of_Fourier_terms_in_PUL_LC_calc=10

  TLM_cell_equivalent_radius_factor   =1.08D0
  Capacitance_equivalent_radius_factor=1.08d0
  Inductance_equivalent_radius_factor =1d0/1.08D0
  
  Max_cable_bundle_diameter_factor    =0.67d0

  CALL write_line('FINISHED: reset_packet_data',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE reset_packet_data
