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
! SUBROUTINE subtract_dc
!
! NAME
!    subtract_dc
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/3/2018 CJS
!
!
SUBROUTINE subtract_dc


USE post_process
USE file_information

IMPLICIT NONE

! local variables
  
  real*8	:: dc
  
  integer	:: time_loop
  
! START

!  write(*,*)'Sum Time Domain Data'
  
  n_functions_of_time=1
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  write(post_process_info_unit,*)'	Time domain functions:'
  
  write(*,*)'File for time domain data:'
  CALL read_time_domain_data(1)
    
  dc=0d0
  
  do time_loop=1,function_of_time(1)%n_timesteps
  
    dc=dc+function_of_time(1)%value(time_loop)
    
  end do ! next timestep
  
  dc=dc/dble(function_of_time(1)%n_timesteps)
  
  do time_loop=1,function_of_time(1)%n_timesteps
  
    function_of_time(1)%value(time_loop)=function_of_time(1)%value(time_loop)-dc
    
  end do ! next time value
  
  CALL write_time_Domain_Data(1)
 
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE subtract_dc
