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
! SUBROUTINE extract_time_domain_data
!
! NAME
!    extract_time_domain_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 17/08/2012 CJS
!
!
SUBROUTINE extract_time_domain_data

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  integer		:: function_number

! START

  n_functions_of_time=1
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  function_number=1
  
  CALL read_Time_Domain_Data(function_number) ! read function of time from output dataset
  
! Write the single data set to file
  
  CALL write_time_domain_data(function_number)
  
  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE extract_time_domain_data
