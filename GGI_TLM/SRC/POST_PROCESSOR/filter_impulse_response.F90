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
! SUBROUTINE filter impulse response
!
! NAME
!    Calculate the impulse response of a filter function given a timestep value
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/2/2018 CJS based on Apply_filter_in_time_domain
!
!
SUBROUTINE filter_impulse_response

USE post_process
USE filter_types		        
USE filter_functions		        
USE file_information

IMPLICIT NONE

! local variables

  integer			:: function_number

! Filter stuff
  character*256			:: filter_file_name
  real*8			:: dt,tmax
  type(Sfilter) 		:: Sfilter1
  type(Zfilter) 		:: Zfilter1
  type(Zfilter_response) 	:: Zfilter_data1

  integer			:: timestep,n_timesteps
  real*8			:: f_input
    
  integer,parameter :: min_timesteps=20

! START

! Allocate and read the time domain input data
  n_functions_of_time=1
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
    
! Allocate and read the filter function to apply

  write(*,*)'Please enter the filename for the filter function'
  read(*,'(A)')filter_file_name
  write(record_user_inputs_unit,'(A)')trim(filter_file_name)

  open(UNIT=local_file_unit, 					&
       FILE=trim(filter_file_name),	    			&
       STATUS='old',						&
       ERR=9000)

  call read_Sfilter(Sfilter1,local_file_unit) ! filter function

  close(UNIT=local_file_unit)
  
  write(post_process_info_unit,*)'	Filter filename:',trim(filter_file_name)

! Calculate the digital filter coefficients using the bilinear transformation
  
  write(*,*)'Enter the timestep for the filter impulse response'
  read(*,*) dt
  write(record_user_inputs_unit,*)dt,'     # timestep for impulse response'
  
  write(*,*)'Enter the maximum time for the filter impulse response'
  read(*,*) tmax
  write(record_user_inputs_unit,*)tmax,'     # maximum time for impulse response'
  
  n_timesteps=NINT(tmax/dt)
  if (n_timesteps.LT.min_timesteps) then
    n_timesteps=min_timesteps
  end if
  
  write(*,*)'Number of timesteps=',n_timesteps
  
  function_of_time(1)%n_timesteps=n_timesteps
  ALLOCATE ( function_of_time(1)%time(1:n_timesteps) )
  ALLOCATE ( function_of_time(1)%value(1:n_timesteps) ) 
  
  Zfilter1=s_to_z(Sfilter1,dt) 

  Zfilter_data1=allocate_Zfilter_response(Zfilter1%a%order,Zfilter1%b%order)
  
! two loops over the impulse response calculation:
! first to calcualte the number of timesteps required and the second to 
! put the impulse response into the function

  
  Zfilter_data1%w(:)=0d0

  do timestep=1,n_timesteps
    
    CALL timeshift_Zfilter(Zfilter_data1)
    
    if (timestep.EQ.1) then
      f_input=1.0
    else
      f_input=0.0
    end if
    
    CALL evaluate_Zfilter(Zfilter1,Zfilter_data1,f_input)
            
    function_of_time(1)%time(timestep)=(timestep-1)*dt
    function_of_time(1)%value(timestep)=Zfilter_data1%f 
     
  end do  ! next timestep
  
! Write the filtered data set to file
  function_number=1
  
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
  
END SUBROUTINE filter_impulse_response
