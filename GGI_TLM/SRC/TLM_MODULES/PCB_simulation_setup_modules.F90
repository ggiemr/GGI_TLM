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

real*8             :: dl     ! cell size for TLM solution

integer,parameter  :: max_gerber_files=20
integer            :: n_gerber_files
character(LEN=256) :: gerber_filename(max_gerber_files)
real*8             :: gerber_z_offset(max_gerber_files)

integer,parameter  :: max_surfaces=100
integer            :: n_surfaces
integer            :: surface_type(max_surfaces)
character(LEN=256) :: surface_filename(max_surfaces)

integer,parameter  :: surface_type_stl=1

integer,parameter  :: max_volumes=100
integer            :: n_volumes
integer            :: volume_type(max_volumes)

integer,parameter  :: volume_type_rectangular_block2=1


END MODULE PCB_simulation
