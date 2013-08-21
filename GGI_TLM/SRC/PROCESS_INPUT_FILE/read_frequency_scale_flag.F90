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
! SUBROUTINE read_frequency_scale_flag
!
! NAME
!     read_frequency_scale_flag
!
! DESCRIPTION
!     read frequency scaling flag and scale factor
!
! Example packet:
!
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 12/02/2013 CJS
!
!
SUBROUTINE read_frequency_scale_flag

USE TLM_general
USE file_information

IMPLICIT NONE

! local variables

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_frequency_scale_flag',0,output_to_screen_flag)

  read(input_file_unit,*,err=9000)frequency_scale_flag
  read(input_file_unit,*,err=9000)frequency_scale

  CALL write_line('FINISHED: read_frequency_scale_flag',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error reading frequency_scale_flag:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_frequency_scale_flag
