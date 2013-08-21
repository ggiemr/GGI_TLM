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
! MODULE FF_parameters
! MODULE FF_general
! MODULE FF_input_data
! MODULE FF_filters
! MODULE FF_file_stuff

!
! NAME
!     MODULE FF_parameters
!
! DESCRIPTION
!     general parameters relating to the vector fitting process
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 10/12/12 CJS based on HIRF-SE and Flaviir work
!
!
MODULE FF_parameters

USE filter_types

IMPLICIT NONE
              
  integer,parameter	:: max_iterations=10
  
  integer,parameter	:: max_opt_iterations=1000
  real*8,parameter 	:: opt_accuracy=1D-8
  
  integer,parameter 	:: dielectric_material =1	    
  integer,parameter	:: magnetic_material   =2	    
  integer,parameter	:: thin_layer	   =3	    
  integer,parameter	:: impedance	   =4

END MODULE FF_parameters

!
! NAME
!     MODULE FF_general
!
! DESCRIPTION
!     general data relating to the filter fitting process
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 7/12/10 CJS based on Flaviir work
!
!
MODULE FF_general

IMPLICIT NONE

! problem_name

  character*256 :: FF_name
  
  integer	:: fit_type
  
  integer	:: order
  real*8	:: fscale
  
  logical	:: stabilise_input_data_flag
  logical	:: stabilise_filter_flag
  logical	:: optimise_filter_flag
  
  logical	:: ff_output_to_screen
  
  logical	:: unstable_input_data
  logical	:: stable_filter
     
END MODULE FF_general

!
! NAME
!     MODULE FF_input_data
!
! DESCRIPTION
!     input frequency domain data 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 7/12/10 CJS based on Flaviir work
!
!
MODULE FF_input_data

IMPLICIT NONE

  integer			:: n_functions
  integer			:: n_values

  real*8,allocatable		:: frequency(:)
  real*8,allocatable		:: normalised_frequency(:)
  complex*16,allocatable	:: s(:)
  complex*16,allocatable 	:: value(:,:)
  
  real*8			:: fnorm
  real*8			:: wnorm
  
  integer,parameter		:: z11=1
  integer,parameter		:: z12=2
  integer,parameter		:: z21=3
  integer,parameter		:: z22=4
   
END MODULE FF_input_data

!
! NAME
!     MODULE FF_testing_data
!
! DESCRIPTION
!     input frequency domain testing data 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 7/12/10 CJS based on Flaviir work
!
!
MODULE FF_testing_data

IMPLICIT NONE

  integer 			:: n_testing_frequencies
  real*8,allocatable		:: testing_frequency(:)
  complex*16,allocatable	:: testing_s(:)
  
END MODULE FF_testing_data
!
! NAME
!     MODULE FF_filters
!
! DESCRIPTION
!     general data relating to the vector fitting process
!	
!     
! COMMENTS
!     
!     uses HIRF-SE filter library
!
!
! HISTORY
!
!     started 7/12/10 CJS
!
!
MODULE FF_filters

USE filter_types

IMPLICIT NONE
  
  integer			:: n_filters
  integer			:: n_sigma
  
  TYPE(Sfilter)   ,allocatable	:: filter_S(:)
  TYPE(Sfilter_PZ),allocatable	:: filter_S_PZ(:)
  TYPE(Sfilter_PR),allocatable	:: filter_S_PR(:)
  real*8,allocatable		:: filter_sigma(:)  
  
  real*8			:: Mean_square_error
  
END MODULE FF_filters

!
! NAME
!     MODULE FF_file_stuff
!
! DESCRIPTION
!     general parameters relating to the vector fitting file information
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 7/12/10 CJS based on Flaviir work
!
!
MODULE FF_file_stuff

IMPLICIT NONE
!
       integer,parameter		:: initial_data_file_unit    =10
       integer,parameter		:: input_data_file_unit      =11
       integer,parameter 		:: filter_fit_data_file_unit =12
       integer,parameter 		:: filter_file_unit          =13
       integer,parameter 		:: testing_frequency_file_unit=14
!
       character(LEN=10),parameter 	:: dielectric_initial_data_extension    ='.eps.fdmat'
       character(LEN=9) ,parameter 	:: magnetic_initial_data_extension      ='.mu.fdmat'
       character(LEN=10),parameter 	:: thin_layer_z11_initial_data_extension='.z11.fdmat'
       character(LEN=10),parameter 	:: thin_layer_z12_initial_data_extension='.z12.fdmat'
       character(LEN=10),parameter 	:: thin_layer_z21_initial_data_extension='.z21.fdmat'
       character(LEN=10),parameter 	:: thin_layer_z22_initial_data_extension='.z22.fdmat'
       character(LEN=8) ,parameter 	:: impedance_initial_data_extension     ='.Z.fdmat'
!       
       character(LEN=13),parameter 	:: dielectric_input_data_extension         ='.eps.fd_input'
       character(LEN=12),parameter 	:: magnetic_input_data_extension           ='.mu.fd_input'
       character(LEN=13),parameter 	:: thin_layer_z11_input_data_extension     ='.z11.fd_input'
       character(LEN=13),parameter 	:: thin_layer_z12_input_data_extension     ='.z12.fd_input'
       character(LEN=13),parameter 	:: thin_layer_z21_input_data_extension     ='.z21.fd_input'
       character(LEN=13),parameter 	:: thin_layer_z22_input_data_extension     ='.z22.fd_input'
       character(LEN=11),parameter 	:: impedance_input_data_extension          ='.Z.fd_input'
!      
       character(LEN=13),parameter 	:: dielectric_trial_extension         ='.eps.fd_trial'
       character(LEN=12),parameter 	:: magnetic_trial_extension           ='.mu.fd_trial'
       character(LEN=13),parameter 	:: thin_layer_z11_trial_extension     ='.z11.fd_trial'
       character(LEN=13),parameter 	:: thin_layer_z12_trial_extension     ='.z12.fd_trial'
       character(LEN=13),parameter 	:: thin_layer_z21_trial_extension     ='.z21.fd_trial'
       character(LEN=13),parameter 	:: thin_layer_z22_trial_extension     ='.z22.fd_trial'
       character(LEN=11),parameter 	:: impedance_trial_extension          ='.Z.fd_trial'
!                            
       character(LEN=5),parameter 	:: dielectric_filter_extension='.vmat'
       character(LEN=5),parameter 	:: magnetic_filter_extension  ='.vmat'
       character(LEN=5),parameter 	:: thin_layer_filter_extension='.smat'
       character(LEN=2),parameter 	:: impedance_filter_extension ='.Z'
       
       character(LEN=26),parameter 	:: testing_frequency_filename ='testing_frequencies.fd_out'

END MODULE FF_file_stuff
