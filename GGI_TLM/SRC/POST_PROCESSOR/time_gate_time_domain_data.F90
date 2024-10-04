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
! SUBROUTINE time_gate_time_domain_data
!
! NAME
!    time_gate_time_domain_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/09/2024 CJS based on extract_time_domain_data.F90
!
!
SUBROUTINE time_gate_time_domain_data

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  integer		:: function_number
  
  real*8	:: tmin,tmax,t
  integer	:: time_loop
  
! START

  n_functions_of_time=1
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  function_number=1
  
  CALL read_Time_Domain_Data(function_number) ! read function of time from output dataset
  
  write(*,*)'Enter minimum time to output, tmin'
  read(*,*)tmin
  write(*,*)'Enter maximum time to output, tmax'
  read(*,*)tmax
  
  write(record_user_inputs_unit,'(E16.7,A)')tmin,' tmin'
  write(record_user_inputs_unit,'(E16.7,A)')tmax,' tmax'

  do time_loop=1,function_of_time(function_number)%n_timesteps
  
    t=function_of_time(function_number)%time(time_loop)
    
    if ( (t.LT.tmin).OR.(t.GT.tmax) ) then
      function_of_time(function_number)%value(time_loop)=0d0
    end if
    
  end do ! next time value
  
! Write the single data set to file
  
  CALL write_time_domain_data(function_number)
  
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE time_gate_time_domain_data
