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
! MODULE geometry_types
! MODULE geometry
!
! NAME
!     MODULE geometry_types
!
! DESCRIPTION
!     data relating to types of geometric entity
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/08/2012 CJS
!

MODULE geometry_types

! Real coordinate types : point, line, triangle, tet

TYPE::xyz

  REAL*8		:: x
  REAL*8		:: y
  REAL*8		:: z
 
END TYPE xyz

TYPE::xyz_line

  type(xyz),dimension(2) 	:: end
 
END TYPE xyz_line

TYPE::xyz_triangle

  type(xyz),dimension(3) 	:: vertex
 
END TYPE xyz_triangle

TYPE::xyz_tet

  type(xyz),dimension(4) 	:: vertex
 
END TYPE xyz_tet

TYPE::ijk

  integer		:: i
  integer		:: j
  integer		:: k
 
END TYPE ijk

TYPE::cell_point

  type(ijk)	:: cell
  integer	:: point
 
END TYPE cell_point

TYPE::cell_segment

  type(cell_point)	:: segment_point(2)
 
END TYPE cell_segment

TYPE::transformation_type

  character(len=256)  	:: trans_type
  INTEGER	    	:: trans_number
  REAL*8		:: parameters(6) 
 
END TYPE transformation_type

TYPE::volume_type

  character(len=256)		:: volume_type_string
  character(len=256)		:: filename
  integer 			:: volume_type
  integer 			:: volume_number
  integer 			:: volume_material_number
  integer 			:: n_volume_parameters
  real*8			:: volume_parameters(20)
  type(transformation_type) 	:: trans
  integer			:: number_of_tets
  type(xyz_tet),allocatable	:: tet_list(:)
  integer			:: number_of_cells
  type(cell_point),allocatable	:: cell_list(:)

END TYPE volume_type

TYPE::surface_type

  character(len=256)		:: surface_type_string
  character(len=256)		:: filename
  integer 			:: surface_type
  integer 			:: surface_number
  integer 			:: surface_material_number
  integer 			:: n_surface_parameters
  real*8			:: surface_parameters(20)
  type(transformation_type) 	:: trans
  integer			:: number_of_triangles
  type(xyz_triangle),allocatable:: triangle_list(:)
  integer			:: number_of_faces
  type(cell_point),allocatable	:: face_list(:)
  
  real*8			:: mesh_xmin,mesh_xmax
  real*8			:: mesh_ymin,mesh_ymax
  real*8			:: mesh_zmin,mesh_zmax
  
  integer			:: mesh_cell_xmin,mesh_cell_xmax
  integer			:: mesh_cell_ymin,mesh_cell_ymax
  integer			:: mesh_cell_zmin,mesh_cell_zmax

END TYPE surface_type

TYPE::line_type

  character(len=256)		:: line_type_string
  character(len=256)		:: filename
  integer 			:: GiD_line_number
  integer 			:: line_type
  integer 			:: line_number
  integer 			:: end_connection(2)
  integer 			:: n_line_parameters
  real*8			:: line_parameters(20)
  type(transformation_type) 	:: trans
  integer			:: number_of_line_segments
  type(xyz_line),allocatable	:: line_segment_list(:)
  integer			:: number_of_cell_segments
  type(cell_segment),allocatable:: cell_segment_list(:)

END TYPE line_type

TYPE::point_type

  type(xyz)			:: point
  type(ijk)			:: cell
  type(cell_point)		:: face
  type(transformation_type) 	:: trans

END TYPE point_type

! VOLUME OBJECT TYPES

  integer,parameter 	  :: volume_type_rectangular_block=1
  integer,parameter 	  :: volume_type_cylinder=2
  integer,parameter 	  :: volume_type_sphere=3
  integer,parameter 	  :: volume_type_rectangular_block2=4
  integer,parameter 	  :: volume_type_tet=5
  integer,parameter 	  :: volume_type_vtk_tet=6
  integer,parameter 	  :: volume_type_pyramid=7
  integer,parameter 	  :: volume_type_pyramid_ram=8
  integer,parameter	  :: volume_type_tet_mesh=9

! SURFACE OBJECT TYPES

  integer,parameter	:: surface_type_rectangular_block=1
  integer,parameter 	:: surface_type_cylinder=2
  integer,parameter	:: surface_type_sphere=3
  integer,parameter	:: surface_type_rectangle=4
  integer,parameter	:: surface_type_circle=5
  integer,parameter	:: surface_type_quad=6
  integer,parameter	:: surface_type_rectangular_block2=7
  integer,parameter	:: surface_type_xplane=8
  integer,parameter	:: surface_type_yplane=9
  integer,parameter	:: surface_type_zplane=10
  integer,parameter	:: surface_type_triangle=11
  integer,parameter	:: surface_type_triangulated_surface=12
  integer,parameter	:: surface_type_vtk_triangulated_surface=13

! LINE OBJECT TYPES

  integer,parameter 	  :: line_type_straight_line=1  
  integer,parameter 	  :: line_type_straight_line2=2 
  integer,parameter 	  :: line_type_arc=3

! LINE END CONNECTION TYPES

  integer,parameter 	  :: line_free=0
  integer,parameter 	  :: line_line=1
  integer,parameter 	  :: line_surface=2

END MODULE geometry_types
!
! NAME
!     MODULE geometry
!
! DESCRIPTION
!     data relating to geometry
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
MODULE geometry

USE geometry_types

IMPLICIT NONE
  
  integer :: n_volumes
  type(volume_type),allocatable 	:: problem_volumes(:)
  
  integer :: n_surfaces
  type(surface_type),allocatable 	:: problem_surfaces(:)
  
  integer :: n_lines
  type(line_type),allocatable 		:: problem_lines(:)
  
  integer :: n_points
  type(point_type),allocatable 		:: problem_points(:)
  
END MODULE geometry
