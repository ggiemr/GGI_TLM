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
! MODULE Mlayer_module
! MODULE Mlayer_filter_module

!
! NAME
!     MODULE Mlayer_module
!
! DESCRIPTION
!     data relating to the M layer S parameter calculation
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
MODULE Mlayer_module

USE filter_types
  
! parameters
  
  integer TE
  integer TM
  
  parameter (TE=1)
  parameter (TM=2)
  
  integer material
  integer thin_layer
  
  parameter (material=1)
  parameter (thin_layer=2)
  
  logical overflow_error


! solution data

  integer m
  
  real*8 angle
  real*8 angle_rad
  
  character*2 polarisation_in
  integer polarisation
  
  real*8 w,f,fmin,fmax,fstep
  
  integer,allocatable :: layer_type(:)
  real*8,allocatable  :: layer_thickness(:)
  
TYPE::material_type

  INTEGER 		:: layer_number
  REAL*8 		:: fmin,fmax
  TYPE(Sfilter)		:: eps_S
  TYPE(Sfilter)		:: mu_S
  REAL*8		:: sigma_e
  REAL*8		:: sigma_m
   
END TYPE material_type

 
TYPE::thin_layer_type

  INTEGER 		:: layer_number
  REAL*8 		:: fmin,fmax
  TYPE(Sfilter)		:: Z11_S,Z12_S,Z21_S,Z22_S
   
END TYPE thin_layer_type

  type(material_type),allocatable  :: material_list(:)
  type(thin_layer_type),allocatable :: thin_layer_list(:)
  complex*16,allocatable	:: ABCD(:,:,:)
 

END MODULE Mlayer_module

!
! NAME
!     MODULE Mlayer_file_module
!
! DESCRIPTION
!     data relating to the M layer S parameter calculation files
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
MODULE Mlayer_file_module

  character*4	:: input_file_extn
  character*5	:: info_file_extn
  character*4	:: S11_file_extn
  character*4	:: S12_file_extn
  character*4	:: S21_file_extn
  character*4	:: S22_file_extn
  character*10	:: Z11_file_extn
  character*10	:: Z12_file_extn
  character*10	:: Z21_file_extn
  character*10	:: Z22_file_extn
  
  character*5	:: power_file_extn
  
  character*4	:: material_file_extn
  character*5	:: thin_layer_file_extn
  
  parameter (input_file_extn=".inp")
  parameter (info_file_extn=".info")
  parameter (S11_file_extn=".S11")
  parameter (S12_file_extn=".S12")
  parameter (S21_file_extn=".S21")
  parameter (S22_file_extn=".S22")
  parameter (Z11_file_extn=".z11.fdmat")
  parameter (Z12_file_extn=".z12.fdmat")
  parameter (Z21_file_extn=".z21.fdmat")
  parameter (Z22_file_extn=".z22.fdmat")
  parameter (power_file_extn=".Pout")
  parameter (material_file_extn=".mat")
  parameter (thin_layer_file_extn=".tmat")
  
  integer	:: input_file_unit
  integer	:: info_file_unit
  integer	:: S11_file_unit
  integer	:: S12_file_unit
  integer	:: S21_file_unit
  integer	:: S22_file_unit
  integer	:: Z11_file_unit
  integer	:: Z12_file_unit
  integer	:: Z21_file_unit
  integer	:: Z22_file_unit
  
  integer	:: power_file_unit
  
  integer	:: material_file_unit
  integer	:: thin_layer_file_unit
  
  parameter (input_file_unit=9)
  parameter (info_file_unit=10)
  parameter (S11_file_unit=11)
  parameter (S12_file_unit=12)
  parameter (S21_file_unit=13)
  parameter (S22_file_unit=14)
  parameter (Z11_file_unit=15)
  parameter (Z12_file_unit=16)
  parameter (Z21_file_unit=17)
  parameter (Z22_file_unit=18)
  parameter (power_file_unit=19)
  parameter (material_file_unit=20)
  parameter (thin_layer_file_unit=21)

  character*256	:: problem_name
  integer	:: problem_name_length    
  
  character*256	:: inp_filename
  character*256	:: info_filename
  character*256	:: S11_filename
  character*256	:: S12_filename
  character*256	:: S21_filename
  character*256	:: S22_filename
  character*256	:: Z11_filename
  character*256	:: Z12_filename
  character*256	:: Z21_filename
  character*256	:: Z22_filename
  
  character*256	:: power_filename
  
  character*256	:: material_filename
  character*256	:: thin_layer_filename

END MODULE Mlayer_file_module
