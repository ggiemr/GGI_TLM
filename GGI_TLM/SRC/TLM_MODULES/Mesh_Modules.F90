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
! MODULE cell_parameters
! MODULE mesh
! MODULE local_grid
!
! NAME
!     MODULE cell_parameters
!
! DESCRIPTION
!     data relating to TLM cells
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/08/2012 CJS
!
MODULE cell_parameters

IMPLICIT NONE

  integer,parameter	:: face_xmin=1
  integer,parameter	:: face_ymin=2
  integer,parameter	:: face_zmin=3
  integer,parameter	:: face_xmax=4
  integer,parameter	:: face_ymax=5
  integer,parameter	:: face_zmax=6
  integer,parameter	:: centre=7
  
  character(len=6),parameter	::face_string(7)=   &
  (/'  xmin','  ymin','  zmin','  xmax','  ymax','  zmax','centre'/)
  
  character(len=6),parameter	:: xmin_string='xmin'
  character(len=6),parameter	:: xmax_string='xmax'
  character(len=6),parameter	:: ymin_string='ymin'
  character(len=6),parameter	:: ymax_string='ymax'
  character(len=6),parameter	:: zmin_string='zmin'
  character(len=6),parameter	:: zmax_string='zmax'
  character(len=6),parameter	:: centre_string='centre'

END MODULE cell_parameters
!
! NAME
!     MODULE mesh
!
! DESCRIPTION
!     data relating to mesh
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
MODULE mesh

USE cell_parameters
USE TLM_general

IMPLICIT NONE

! TLM port numbering

  integer,parameter :: Vy_xmin=1
  integer,parameter :: Vz_xmin=2
  integer,parameter :: Vy_xmax=3
  integer,parameter :: Vz_xmax=4
  integer,parameter :: Vx_ymin=5
  integer,parameter :: Vz_ymin=6
  integer,parameter :: Vx_ymax=7
  integer,parameter :: Vz_ymax=8
  integer,parameter :: Vx_zmin=9
  integer,parameter :: Vy_zmin=10
  integer,parameter :: Vx_zmax=11
  integer,parameter :: Vy_zmax=12

! TLM voltage pulse array

  real*8,allocatable  :: V(:,:,:,:)   ! V(port,i,j,k)
  
! TLM voltages for parallel transfer
  real*8,allocatable  :: Vi_zmin(:)
  real*8,allocatable  :: Vi_zmax(:)
  real*8,allocatable  :: Vr_zmin(:)
  real*8,allocatable  :: Vr_zmax(:)  
  
! cell centre update array: indicates
! what is to be done at a particular cell

  integer		:: number_of_cell_centre_codes
  integer,allocatable	:: cell_centre_update_code(:)   
  
! internal boundary condition array: indicates
! what is to be done at a particular internal boundary face

  integer		:: number_of_face_codes
  integer,allocatable	:: face_update_code(:)   


! temporary mesh arrays to hold material codes  
  integer,allocatable	:: local_surface_material(:,:,:,:)
  integer,allocatable	:: local_cell_material(:,:,:)
  
! temporary mesh arrays to hold output codes  
  integer,allocatable	:: local_surface_output(:,:,:,:)
  integer,allocatable	:: local_cell_output(:,:,:)
  
! temporary mesh arrays to hold cable codes  
  integer,allocatable	:: local_surface_cable(:,:,:,:)
  integer,allocatable	:: local_cell_cable(:,:,:)
  
! temporary mesh arrays to hold excitation codes  
  integer,allocatable	:: local_surface_excitation(:,:,:,:)
  integer,allocatable	:: local_cell_excitation(:,:,:)
  
! code lookup tables  
  integer		:: n_special_cells
  integer,allocatable	:: cell_update_code_to_material_data(:,:)
  integer,allocatable	:: cell_update_code_to_cable_cell_number(:)
  integer,allocatable	:: cell_update_code_to_excitation_number(:)
  integer,allocatable	:: cell_update_code_to_output_number(:)

  integer		:: n_special_faces
  integer,allocatable	:: face_update_code_to_material_data(:,:)
  integer,allocatable	:: face_update_code_to_cable_number(:)
  integer,allocatable	:: face_update_code_to_excitation_number(:)
  integer,allocatable	:: face_update_code_to_output_number(:)

END MODULE mesh
!
! NAME
!     MODULE local_grid
!
! DESCRIPTION
!     data relating to the local grid used for meshing trinagles and tets
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/05/2013 CJS
!
MODULE local_mesh

USE cell_parameters
USE TLM_general

IMPLICIT NONE

! numbering

  integer,parameter :: local_corner=1
  integer,parameter :: local_xface =2
  integer,parameter :: local_yface =3
  integer,parameter :: local_zface =4
  integer,parameter :: local_centre=5
  integer,parameter :: local_xedge =6
  integer,parameter :: local_yedge =7
  integer,parameter :: local_zedge =8
  
  integer,parameter :: not_set=-99999
  
  real*8 		:: local_xmin,local_ymin,local_zmin
  
  integer 		:: local_ixmin,local_iymin,local_izmin
  integer 		:: local_ixmax,local_iymax,local_izmax
  
  integer,allocatable	:: local_grid(:,:,:,:)
  
  integer,allocatable	:: local_grid_tet(:,:,:,:)
  
  integer		:: tot_n_faces
  integer		:: n_faces_set
  
  integer		:: n_x_faces
  integer,allocatable	:: local_x_faces(:,:)
  
  integer		:: n_y_faces
  integer,allocatable	:: local_y_faces(:,:)
  
  integer		:: n_z_faces
  integer,allocatable	:: local_z_faces(:,:)

END MODULE local_mesh
