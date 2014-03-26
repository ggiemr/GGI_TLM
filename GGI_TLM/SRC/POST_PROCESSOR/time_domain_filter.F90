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
! SUBROUTINE Apply_filter_in_time_domain
!
! NAME
!    Apply_filter_in_time_domain
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/1/2013 CJS
!
!
SUBROUTINE Apply_filter_in_time_domain

USE post_process
USE filter_types		        
USE filter_functions		        
USE file_information

IMPLICIT NONE

! local variables

  integer			:: function_number

! Filter stuff
  character*256			:: filter_file_name
  real*8			:: dt
  type(Sfilter) 		:: Sfilter1
  type(Zfilter) 		:: Zfilter1
  type(Zfilter_response) 	:: Zfilter_data1

  integer			:: timestep,n_timesteps
  real*8			:: f_input

! START

! Allocate and read the time domain input data
  n_functions_of_time=2
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  function_number=1
  
  CALL read_Time_Domain_Data(function_number) ! read function of time from output dataset
  
! Allocate and read the filter function to apply

  read(*,'(A)')filter_file_name
  write(record_user_inputs_unit,'(A)')trim(filter_file_name)

  open(UNIT=local_file_unit, 					&
       FILE=trim(filter_file_name),	    			&
       STATUS='old',						&
       ERR=9000)

  call read_Sfilter(Sfilter1,local_file_unit) ! filter function

  close(UNIT=local_file_unit)

! Calculate the digital filter coefficients using the bilinear transformation
  
  dt=function_of_time(1)%time(2)-function_of_time(1)%time(1)
  Zfilter1=s_to_z(Sfilter1,dt) 

  Zfilter_data1=allocate_Zfilter_response(Sfilter1%a%order,Sfilter1%b%order)

! allocate memory for the time domain response function

  n_timesteps=function_of_time(1)%n_timesteps
  
  write(*,*)'Allocating time response',n_timesteps
  
  function_of_time(2)%n_timesteps=n_timesteps
  ALLOCATE ( function_of_time(2)%time(1:n_timesteps) )
  ALLOCATE ( function_of_time(2)%value(1:n_timesteps) )
  
! Apply the digital filter to the time domain data
  
  Zfilter_data1%w(:)=0d0
  
  do timestep=1,n_timesteps
  
    CALL timeshift_Zfilter(Zfilter_data1)
    
    f_input=function_of_time(1)%value(timestep)
    
    CALL evaluate_Zfilter(Zfilter1,Zfilter_data1,f_input)
    
    function_of_time(2)%time(timestep)=function_of_time(1)%time(timestep)
    function_of_time(2)%value(timestep)=Zfilter_data1%f  
    
  end do ! next sample
  
! Write the filtered data set to file
  function_number=2
  
  CALL write_time_domain_data(function_number)

! deallocate memory  
  CALL Deallocate_post_data()
  CALL deallocate_Sfilter( Sfilter1 )
  CALL deallocate_Zfilter( Zfilter1 )
  CALL deallocate_Zfilter_data( Zfilter_data1 )
  
  RETURN
  
9000 CALL write_line('Error reading filter function',0,.TRUE.)
     CALL write_line('Problem opening filter file:'//trim(filter_file_name),0,.TRUE.)
     STOP
  
END SUBROUTINE Apply_filter_in_time_domain
