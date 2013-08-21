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
! MODULE Cables
!
! NAME
!     MODULE Cables
!
! DESCRIPTION
!     general data relating to individual cables
!     and data relating to cable bundles resulting from multiple cables in a cell
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!
!
MODULE Cables

USE Geometry_types
USE Filter_types

IMPLICIT NONE

TYPE::cable_geometry_type 

  character*256		:: cable_geometry_type_string  
  integer		:: cable_geometry_type
  integer		:: n_parameters
  real*8,allocatable	:: parameters(:)
  integer		:: n_conductors
  integer		:: n_external_conductors
  real*8,allocatable	:: external_conductor_xc(:)
  real*8,allocatable	:: external_conductor_yc(:)
  real*8,allocatable	:: external_conductor_radius(:)
  real*8,allocatable	:: external_dielectric_radius(:)
  real*8,allocatable	:: external_dielectric_permittivity(:)
  integer		:: n_shielded_conductors
  real*8,allocatable	:: shielded_conductor_xc(:)
  real*8,allocatable	:: shielded_conductor_yc(:)
  real*8,allocatable	:: shielded_conductor_radius(:)
  real*8,allocatable	:: shielded_dielectric_radius(:)
  real*8,allocatable	:: shielded_dielectric_permittivity(:)
  real*8		:: cable_offset_radius
  integer,allocatable	:: SC(:)
  integer,allocatable	:: Tv(:,:)
  integer,allocatable	:: Ti(:,:)
  real*8,allocatable	:: L_internal(:,:)
  real*8,allocatable	:: C_internal(:,:)
  real*8,allocatable	:: R_internal(:,:)
  
  integer			:: n_filters
  integer,allocatable		:: filter_number(:,:)  
  type(Sfilter),allocatable	:: Sfilter(:)
  type(Zfilter),allocatable	:: Zfilter(:)
  real*8,allocatable		:: Z_f(:)
  
END TYPE cable_geometry_type

  integer,parameter	:: cable_geometry_type_cylindrical	=1
  integer,parameter	:: cable_geometry_type_FD_cylindrical	=2
  integer,parameter	:: cable_geometry_type_coaxial		=3
  integer,parameter	:: cable_geometry_type_FD_coaxial	=4
  integer,parameter	:: cable_geometry_type_ribbon		=5
  integer,parameter	:: cable_geometry_type_FD_ribbon	=6

  integer			:: n_cable_geometries
  type(cable_geometry_type),allocatable	::cable_geometry_list(:)   
  
TYPE::cable_type  

  integer			:: cable_geometry_number
  integer			:: n_conductors
  integer			:: n_lines
  integer,allocatable		:: line_list(:)
  integer			:: junction_1
  integer			:: junction_2
  integer			:: number_of_cable_segments
  type(cell_segment),allocatable:: cable_segment_list(:)
  integer,allocatable		:: direction_sign_list(:)
  integer,allocatable		:: bundle_segment_list(:)

  integer			:: bundle_junction_1
  integer			:: first_internal_connection_node_in_bundle_junction_1
  integer			:: first_external_conductor_in_bundle_junction_1
  
  integer			:: bundle_junction_2
  integer			:: first_internal_connection_node_in_bundle_junction_2
  integer			:: first_external_conductor_in_bundle_junction_2
   
END TYPE cable_type

  integer			:: n_cables
  type(cable_type),allocatable	:: cable_list(:)
 
TYPE::cable_output_type  

  integer		:: cable_number
  integer		:: closest_point_number
  integer		:: bundle_segment_number
  type(cell_point)	:: output_point
  integer		:: rank
  integer		:: n_conductors
  integer,allocatable	:: conductor_list(:)
   
END TYPE cable_output_type

  integer				:: n_cable_outputs
  type(cable_output_type),allocatable	:: cable_output_list(:)
  
TYPE::P_matrix

  integer,allocatable	::P(:,:)

END TYPE  P_matrix
  
TYPE::integer_vector

  integer,allocatable	::value(:)

END TYPE  integer_vector
  
TYPE::real8_vector

  real*8,allocatable	::value(:)

END TYPE  real8_vector
  
TYPE::cable_junction_type  

  integer				:: point_number
  type(cell_point)			:: cell_point
  integer				:: junction_type
  integer				:: n_internal_connection_nodes
  integer				:: number_of_cables
  integer,allocatable			:: cable_list(:)
  integer,allocatable			:: cable_end_list(:)
  integer,allocatable			:: n_external_conductors(:)
  type(P_matrix),allocatable		:: Pmatrix(:)
  type(integer_vector),allocatable	:: excitation_function(:)
  type(real8_vector),allocatable	:: resistance(:)
  integer,allocatable			:: BC(:)
  
  integer				:: number_of_internal_impedances
  integer,allocatable			:: node_1(:)
  integer,allocatable			:: node_2(:)
  type(Sfilter),allocatable		:: Sfilter(:)
  
  integer				:: bundle_junction_number
  integer				:: first_internal_connection_node_in_bundle_junction
  integer				:: first_external_conductor_in_bundle_junction
   
END TYPE cable_junction_type

  integer,parameter			:: junction_type_cell=1
  integer,parameter			:: junction_type_face=2

  integer				:: n_cable_junctions
  type(cable_junction_type),allocatable	:: cable_junction_list(:)
   
  integer			:: total_number_cable_cells
  integer			:: total_number_cable_faces

! CABLE BUNDLE STRUCTURES FOLLOW
  
TYPE::bundle_segment_type  

  type(cell_segment)		:: cable_segment
  integer			:: n_cables
  integer,allocatable		:: cable_list(:)
  integer,allocatable		:: direction_sign_list(:)
  integer			:: n_conductors
  integer			:: bundle_segment_geometry
  
  real*8,allocatable		:: xc(:)
  real*8,allocatable		:: yc(:)
  real*8,allocatable		:: rc(:)
  real*8,allocatable		:: ri(:)
  
  real*8,allocatable		:: L(:,:)
  real*8,allocatable		:: C(:,:)
  real*8,allocatable		:: R(:,:)
  integer,allocatable		:: Tv(:,:)
  integer,allocatable		:: Ti(:,:)
  integer,allocatable		:: SC(:)
  real*8,allocatable		:: Zlink(:,:)
  real*8,allocatable		:: Ylink(:,:)
  real*8,allocatable		:: ZLstub(:,:)
  real*8,allocatable		:: Yf(:,:)
  real*8,allocatable		:: Vlink(:)
  real*8,allocatable		:: VLstub(:)
  real*8,allocatable		:: Vsource(:)
  real*8,allocatable		:: Iw_centre(:)
  real*8,allocatable		:: Iw_face(:)
  integer,allocatable		:: excitation_function(:)
  
  real*8			:: cable_bundle_radius
  real*8			:: TLM_cell_equivalent_radius
  real*8			:: TLM_reference_radius_rL
  real*8			:: TLM_reference_radius_rC
  
  integer				:: n_filters
  integer,allocatable			:: filter_number(:,:)  
  type(Sfilter),allocatable		:: Sfilter(:)
  TYPE(Zfilter),allocatable		:: Zfilter(:)
  real*8,allocatable			:: Z_f(:)
  integer				:: n_filter_data
  type(Zfilter_response),allocatable	:: Zfilter_data(:) 
   
END TYPE bundle_segment_type

  integer				:: n_bundle_segments
  type(bundle_segment_type),allocatable	:: bundle_segment_list(:)
  
TYPE::bundle_segment_geometry_type  

  integer			:: n_cables
  integer,allocatable		:: cable_list(:)
  integer			:: n_conductors
  
  real*8,allocatable		:: xc(:)
  real*8,allocatable		:: yc(:)
  real*8,allocatable		:: rc(:)
  real*8,allocatable		:: ri(:)
  
  real*8			:: cable_bundle_radius
  real*8			:: TLM_reference_radius_rL
  real*8			:: TLM_reference_radius_rC
  
  real*8,allocatable		:: L(:,:)
  real*8,allocatable		:: C(:,:)
  real*8,allocatable		:: R(:,:)
  integer,allocatable		:: Tv(:,:)
  integer,allocatable		:: Ti(:,:)
  integer,allocatable		:: SC(:)
  real*8,allocatable		:: Zlink(:,:)
  real*8,allocatable		:: Ylink(:,:)
  real*8,allocatable		:: ZLstub(:,:)
  real*8,allocatable		:: Yf(:,:)
  
  integer				:: n_filters
  integer,allocatable			:: filter_number(:,:)  
  type(Sfilter),allocatable		:: Sfilter(:)
  TYPE(Zfilter),allocatable		:: Zfilter(:)
  real*8,allocatable			:: Z_f(:)
   
END TYPE bundle_segment_geometry_type

  integer				:: n_bundle_segment_geometries
  type(bundle_segment_type),allocatable	:: bundle_segment_geometry_list(:)
 
TYPE::bundle_junction_type 

  type(cell_point)		:: cell_point
  integer			:: n_internal_connection_nodes
  integer			:: internal_connection_node_count
  integer			:: n_segments
  integer,allocatable		:: segment_list(:)
  integer,allocatable		:: n_external_conductors(:)
  integer,allocatable		:: external_conductor_count(:)
  integer,allocatable		:: BC(:)
  type(P_matrix),allocatable	:: P_matrix_list(:)
  
  integer			:: number_of_cable_junctions
  integer,allocatable 		:: cable_junction_list(:)
  
  real*8,allocatable		:: Yf(:,:)
  
  integer			:: n_internal_impedance_filters
  type(Sfilter),allocatable	:: Sfilter(:)
  TYPE(Zfilter),allocatable	:: Zfilter(:)
  real*8,allocatable		:: Z_f(:)
  type(Zfilter_response),allocatable	:: Zfilter_data(:) 
   
END TYPE bundle_junction_type

  integer					:: n_cell_centre_junctions
  type(bundle_junction_type),allocatable	:: cell_centre_junction_list(:)

  integer					:: n_face_junctions
  type(bundle_junction_type),allocatable	:: face_junction_list(:)

! Parallel data passing variables       
         
  integer  		:: n_zmin_segments_send
  integer  		:: n_zmin_reals_send
  integer,allocatable 	:: zmin_segment_list_send(:)
  
  integer  		:: n_zmin_segments_rcv
  integer  		:: n_zmin_reals_rcv
  integer,allocatable 	:: zmin_segment_list_rcv(:)
  
  integer  		:: n_zmax_segments_send
  integer  		:: n_zmax_reals_send
  integer,allocatable 	:: zmax_segment_list_send(:)
  
  integer  		:: n_zmax_segments_rcv
  integer  		:: n_zmax_reals_rcv
  integer,allocatable 	:: zmax_segment_list_rcv(:)
  
  real*8,allocatable 	:: wire_Vi_zmin(:)
  real*8,allocatable 	:: wire_Vi_zmax(:)
  real*8,allocatable 	:: wire_Vr_zmin(:)
  real*8,allocatable 	:: wire_Vr_zmax(:)
       
! Local cable stuff

END MODULE Cables
