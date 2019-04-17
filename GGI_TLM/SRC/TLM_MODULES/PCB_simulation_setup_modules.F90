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
!     MODULE PCB_simulation
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 14/3/19 CJS
!     
!
MODULE PCB_simulation

IMPLICIT NONE

SAVE

! Mesh stuff

real*8             :: dl     ! cell size for TLM solution

real*8             :: xmin,xmax,ymin,ymax,zmin,zmax     ! TLM problem space dimensions

integer            :: mesh_nx,mesh_ny,mesh_nz                          ! number of cells in each direction

! Mesh

integer,allocatable :: material_mesh(:,:,:,:)

integer,parameter :: centre=1
integer,parameter :: face_xmin=2
integer,parameter :: face_xmax=3
integer,parameter :: face_ymin=4
integer,parameter :: face_ymax=5
integer,parameter :: face_zmin=6
integer,parameter :: face_zmax=7

! Gerber file stuff

integer,parameter  :: max_gerber_files=20
integer            :: n_gerber_files
character(LEN=256) :: gerber_filename(max_gerber_files)
real*8             :: gerber_z_offset(max_gerber_files)
character(LEN=256) :: executable_name
character(LEN=256) :: executable_path

! Surface geometry stuff

integer,parameter  :: max_surfaces=100
integer            :: n_surfaces
integer            :: surface_type(max_surfaces)
character(LEN=256) :: surface_filename(max_surfaces)
real*8             :: surface_z_offset(max_surfaces)
real*8             :: surface_parameters(max_surfaces,1:6)

integer,parameter  :: surface_type_stl=1
integer,parameter  :: surface_type_xplane=2
integer,parameter  :: surface_type_yplane=3
integer,parameter  :: surface_type_zplane=4

! Surface material stuff

integer,parameter  :: max_surface_materials=100
integer            :: n_surface_materials
integer            :: surface_material_type(max_surface_materials)

integer            :: surface_material_to_surface_list(max_surface_materials)         

integer,parameter  :: surface_material_type_PEC=1
integer,parameter  :: surface_material_type_SPICE=2

integer,parameter  :: max_PEC_surfaces=100
integer            :: n_PEC_surfaces
integer            :: PEC_surface_list(max_PEC_surfaces)
integer            :: PEC_surface_orientation_list(max_PEC_surfaces)

integer            :: SPICE_port_list(max_surface_materials)         
integer            :: SPICE_node_list(max_surface_materials,2)         
character(LEN=2)   :: SPICE_port_direction_list(max_surface_materials)

Character(LEN=34)  :: terminal_connection_geometry_filename='component_terminal_connections.stl'

! Volume geometry stuff

integer,parameter  :: max_volumes=100
integer            :: n_volumes
integer            :: volume_type(max_volumes)
integer,parameter  :: volume_type_rectangular_block2=1
real*8             :: volume_parameters(max_volumes,1:6)


! Volume material stuff

integer             :: n_volume_materials

integer             :: volume_material_type(max_volumes)
integer,parameter   :: volume_material_type_DISPERSIVE=1

character(LEN=256)  :: volume_material_name(max_volumes)
integer             :: volume_material_to_volume_list(max_volumes)         

! Dielectric stuff

integer :: n_dielectrics

! Via stuff

integer :: n_vias

! lumped component stuff

integer,parameter :: max_lumped_components=100
integer,parameter :: max_ngspice_nodes=2
integer,parameter :: max_ngspice_terminals=3  ! number of nodes+1
integer,parameter :: max_ngspice_ports=200    ! twice the number of lumped components

integer,parameter :: component_type_one_port_model=1
integer,parameter :: component_type_two_port_model=2

integer :: n_lumped_components
integer :: lumped_component_type(max_lumped_components)

integer :: ngspice_n_ports(max_lumped_components)

integer :: ngspice_n_terminals(max_lumped_components)
real*8  :: ngspice_terminal_list(max_lumped_components,max_ngspice_terminals,3)

integer :: tot_n_ngspice_nodes
integer :: ngspice_n_nodes(max_lumped_components)
integer :: ngspice_node_list(max_lumped_components,max_ngspice_nodes,2)
integer :: ngspice_port_list(max_lumped_components,max_ngspice_nodes)

real*8  :: z_position(max_lumped_components)

integer :: package_type(max_lumped_components)

! port information
integer :: tot_n_ngspice_ports
integer :: ngspice_port_to_node_list(max_ngspice_ports,2)

integer,parameter :: package_type_none=0
integer,parameter :: package_type_rectangular_block=1
integer,parameter :: package_type_cylinder=2

real*8  :: package_parameters(max_lumped_components,1:6)
real*8  :: package_orientation(max_lumped_components)

! additional component stuff

integer :: n_additional_components



END MODULE PCB_simulation
