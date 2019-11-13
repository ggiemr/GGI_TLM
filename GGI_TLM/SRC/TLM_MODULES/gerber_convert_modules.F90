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
!     MODULE gerber_parameters
!     MODULE gerber
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 10/12/19 CJS
!     
!
MODULE gerber_parameters

IMPLICIT NONE

integer,parameter :: line_length=1000
integer,parameter :: maxDcodes=100
integer,parameter :: maxAparams=10                ! note, this includes aperture macros
integer,parameter :: maxAM=10                     ! maximum number of aperture macros
integer,parameter :: maxAM_primitives=100         ! maximum number aperture macro primitives
integer,parameter :: maxAM_modifiers=100          ! maximum number of modifiers in an aperture macro primitive

integer,parameter :: linear=1
integer,parameter :: clockwise=2
integer,parameter :: counterclockwise=3

integer,parameter :: single=1
integer,parameter :: multi=2

integer,parameter :: dark=1
integer,parameter :: clear=2

END MODULE gerber_parameters
!
! _______________________________________________________
!
!
MODULE gerber

USE gerber_parameters

IMPLICIT NONE

SAVE

integer   :: nDcodes
integer   :: Dcode(1:maxDcodes)
character :: Atype(1:maxDcodes)
real*8    :: Aparams(1:maxDcodes,1:maxAparams)
real*8    :: Asize(1:maxDcodes)
integer   :: A_AMnumber(1:maxDcodes)
integer   :: A_nvars(1:maxDcodes)
character*256    :: Avars(1:maxDcodes,1:maxAparams)

integer :: aperture

logical :: region
logical :: first_region_point_set
real*8  :: frx,fry

integer :: operation
integer :: interpolation_mode

integer :: quadrant_mode

integer :: polarity

real*8 :: xmin,xmax,ymin,ymax,dframe,max_Asize
integer :: nx,ny
real*8 :: dl
integer,allocatable :: p(:,:)
integer,allocatable :: reg(:,:)
integer ipx,ipy

integer,allocatable :: ap(:,:)
integer :: anx,any

real*8 :: s           ! scale factor to metres

! Triangulation data
integer :: n_nodes
integer :: n_triangles
real*8,allocatable  :: node_list(:,:)
integer,allocatable :: triangle_to_node_list(:,:)

! Aperture macro data

integer       :: n_AM
character*256 :: AMname(1:maxAM)
integer       :: total_AM_primitives
integer       :: AM_n_primitives(1:maxAM)
integer       :: AM_to_primitive_list(1:maxAM,1:maxAM_primitives)
integer       :: AM_primitive_number(1:maxAM_primitives)
integer       :: n_AM_modifiers(1:maxAM_primitives)
character*256 :: AM_modifiers(1:maxAM_primitives,1:maxAM_modifiers)


logical :: verbose=.TRUE.

logical :: AMflag=.FALSE.
logical :: AMflag_outline=.FALSE.
logical :: AMflag_Moire=.FALSE.
logical :: AMflag_thermal=.FALSE.
logical :: ABflag=.FALSE.
logical :: LMflag=.FALSE.
logical :: LRflag=.FALSE.
logical :: LSflag=.FALSE.
logical :: SRflag=.FALSE.
logical :: TFflag=.FALSE.
logical :: TAflag=.FALSE.
logical :: TOflag=.FALSE.
logical :: TDflag=.FALSE.

END MODULE gerber
