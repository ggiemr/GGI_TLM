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
! MODULE TLM_periodic
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
  
  character(len=256)	:: GGI_TLM_version
  
  character(len=256)	:: GGI_TLM_date
    
  character(len=256)	:: GGI_TLM_compilation_date
 
  logical	:: bicubic_warp_flag
  logical	:: frequency_scale_flag
  real*8	:: frequency_scale

  logical	:: write_info_file=.FALSE.
  
  logical	:: output_to_screen_flag=.TRUE.
  logical	:: timestepping_output_to_screen_flag=.FALSE.
  
  logical	:: write_geometry_vtk_files
  
  integer	:: number_of_warnings
  
  logical	:: periodic_boundary

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
  
  real*8	:: delay_x
  real*8	:: delay_y
  real*8	:: delay_z
  
  integer	:: wrap_array_point
  
  real*8,allocatable	:: Vy_xmin_save(:,:)
  real*8,allocatable	:: Vz_xmin_save(:,:)
  real*8,allocatable	:: Vy_xmax_save(:,:)
  real*8,allocatable	:: Vz_xmax_save(:,:)
  
  real*8,allocatable	:: Vx_ymin_save(:,:)
  real*8,allocatable	:: Vz_ymin_save(:,:)
  real*8,allocatable	:: Vx_ymax_save(:,:)
  real*8,allocatable	:: Vz_ymax_save(:,:)
  
  real*8,allocatable	:: Vx_zmin_save(:,:)
  real*8,allocatable	:: Vy_zmin_save(:,:)
  real*8,allocatable	:: Vx_zmax_save(:,:)
  real*8,allocatable	:: Vy_zmax_save(:,:)
  
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
  
  logical	:: set_random_number_seed
  
  
END MODULE TLM_general
!
! NAME
!     MODULE TLM_periodic
!
! DESCRIPTION
!      Data required for modelling periodic structures. See run_TLM_periodic_BC.F90
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/03/2014 CJS
!
MODULE TLM_periodic

IMPLICIT NONE

  
! periodic structure parameters
  
  real*8  t_cycle
  integer n_cycles
  real*8  tmax_cycle
  integer ntmax_cycle
  integer nt_cycle,ntmin,ntmax,nt_cycle_start
  integer nt_run(4),n_pbc_timesteps(4),n_save(4)
  integer ntdx,ntdy

  integer pbc_xface,pbc_yface
  
  integer psox,psoy,psoz
  
  integer n_xbc_point,n_ybc_point
  
  integer xbc_write_point_ymin,xbc_write_point_ymax
  integer ybc_write_point_xmin,ybc_write_point_xmax
  
  integer xbc_write_point_save_ymin,xbc_write_point_save_ymax 
  integer ybc_write_point_save_xmin,ybc_write_point_save_xmax

! perodic boundary transfer data  

  real*8,allocatable  :: V_pbc_save(:,:,:,:)
  real*8,allocatable  :: V_pbc_x(:,:,:,:)
  real*8,allocatable  :: V_pbc_y(:,:,:,:)
  real*8,allocatable  :: V_pbc_x_save(:,:,:,:)
  real*8,allocatable  :: V_pbc_y_save(:,:,:,:)

  
END MODULE TLM_periodic
