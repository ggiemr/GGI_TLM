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
! SUBROUTINE read_small_to_zero
!
! NAME
!     read_small_to_zero
!
! DESCRIPTION
!     set a flag which sets small numbers to zero in the TLM mesh to avoid problems with denormal numbers
!
! Example packet:
!
!set_small_to_zero
!50                     ! set small numbers to zero every nth timestep
!1E-30                  ! value of small
!
! COMMENTS
!     
!
! HISTORY
!
!     started 16/11/2020 CJS try to eliminate a problem which may be due to denormal numbers
!
!
SUBROUTINE read_small_to_zero

USE TLM_general
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: read_small_to_zero',0,output_to_screen_flag)
  
  set_small_to_zero=.TRUE.

  read(input_file_unit,*,err=9000,end=9000)n_small_to_zero
  
  read(input_file_unit,*,err=9010,end=9010)small_to_zero_value

  CALL write_line('FINISHED: read_small_to_zero',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error reading n_small_to_zero:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9010 CALL write_line('Error reading small_to_zero_value:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
    
END SUBROUTINE read_small_to_zero
