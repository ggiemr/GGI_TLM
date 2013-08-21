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
! MODULE TLM_general
!
! NAME
!     MODULE TLM_general
!
! DESCRIPTION
!     general data relating to the overall problem
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!     Parallel stuff started 21/11/2012 CJS
!
MODULE TLM_general

IMPLICIT NONE

#if defined(MPI)

  INCLUDE "mpif.h"
  
  integer,DIMENSION(MPI_STATUS_SIZE) :: status
  
#endif

  character(len=256)	:: problem_name
  
  character(len=256)	:: solver_version
  
  character(len=256)	:: solver_date
  
  logical	:: bicubic_warp_flag
  logical	:: frequency_scale_flag
  real*8	:: frequency_scale
  
  logical	:: output_to_screen_flag=.TRUE.
  logical	:: timestepping_output_to_screen_flag=.FALSE.
  
  integer	:: number_of_warnings

! number of cells in x, y and z directions  
  integer 	:: nx,ny,nz
  
! parallel stuff
  integer,save 	:: rank=0
  integer,save 	:: np=1
  integer 	:: ierr,ierror

  integer	:: nz1,nz2
  integer	:: nzmin,nzmax

  integer,allocatable	:: cell_rank(:)   
  integer,allocatable	:: cell_face_rank(:,:)   
  
! cell dimension (m)  
  real*8	:: dl

! mesh limits in x, y and z directions  
  real*8 mesh_xmin,mesh_xmax
  real*8 mesh_ymin,mesh_ymax
  real*8 mesh_zmin,mesh_zmax

! outer boundary voltage reflection coefficients  
  real*8	:: R_xmin,R_xmax
  real*8	:: R_ymin,R_ymax
  real*8	:: R_zmin,R_zmax
  
! wrapping boundary conditions on outer boundary

  logical 	:: wrap_x
  logical 	:: wrap_y
  logical 	:: wrap_z
  
  real*8,allocatable	:: Vx_wrap_zmin_send(:,:)
  real*8,allocatable	:: Vy_wrap_zmin_send(:,:)
  real*8,allocatable	:: Vx_wrap_zmax_send(:,:)
  real*8,allocatable	:: Vy_wrap_zmax_send(:,:)
  
  real*8,allocatable	:: Vx_wrap_zmin_rcv(:,:)
  real*8,allocatable	:: Vy_wrap_zmin_rcv(:,:)
  real*8,allocatable	:: Vx_wrap_zmax_rcv(:,:)
  real*8,allocatable	:: Vy_wrap_zmax_rcv(:,:)
  
! general solver parameters 
  real*8	:: simulation_time
  real*8	:: dt
  integer	:: n_timesteps
  
  character(len=2),parameter	:: Ex_string='ex'
  character(len=2),parameter	:: Ey_string='ey'
  character(len=2),parameter	:: Ez_string='ez'
  character(len=2),parameter	:: Hx_string='hx'
  character(len=2),parameter	:: Hy_string='hy'
  character(len=2),parameter	:: Hz_string='hz'
  
  character(len=2),parameter	:: Jx_string='jx'
  character(len=2),parameter	:: Jy_string='jy'
  character(len=2),parameter	:: Jz_string='jz'
  character(len=2),parameter	:: Emagnitude_string='em'
  character(len=2),parameter	:: Hmagnitude_string='hm'
  character(len=2),parameter	:: Jmagnitude_string='jm'
  character(len=2),parameter	:: Power_string='po'

  integer,parameter	:: Ex=1  
  integer,parameter	:: Ey=2  
  integer,parameter	:: Ez=3  
  integer,parameter	:: Hx=4  
  integer,parameter	:: Hy=5  
  integer,parameter	:: Hz=6  

  integer,parameter	:: Jx=7
  integer,parameter	:: Jy=8  
  integer,parameter	:: Jz=9  
  integer,parameter	:: Mx=10  
  integer,parameter	:: My=11 
  integer,parameter	:: Mz=12  
  
  integer,parameter	:: Emagnitude=13
  integer,parameter	:: Hmagnitude=14
  integer,parameter	:: Jmagnitude=15
  integer,parameter	:: Mmagnitude=16
  integer,parameter	:: Power=17

  integer	:: timestep
  
  real*8	:: time
  
! Some solution parameters which may be re-set from the input file

  logical	:: new_mesh_generation

  integer		:: Cable_LC_Correction_type
  integer,parameter 	:: LC_correction_type_geometry_scale=1
  integer,parameter 	:: LC_correction_type_subtract_cell_inductance=2
  
  integer	:: number_of_Fourier_terms_in_PUL_LC_calc

  real*8	:: TLM_cell_equivalent_radius_factor
  real*8	:: Capacitance_equivalent_radius_factor
  real*8	:: Inductance_equivalent_radius_factor
  
  real*8	:: Max_cable_bundle_diameter_factor
  
END MODULE TLM_general
