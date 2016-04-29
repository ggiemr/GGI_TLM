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
! SUBROUTINE read_reduced_c_factor
!
! NAME
!     read_reduced_c_factor
!
! DESCRIPTION
!     read the reduced c factor to apply
!
! Example packet:
!
!reduced_c_factor
!10                  (real)
!
! COMMENTS
!     
!
! HISTORY
!
!     started 7/04/2016 CJS
!
!
SUBROUTINE read_reduced_c_factor

USE TLM_general
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: read_reduced_c_factor',0,output_to_screen_flag)

  read(input_file_unit,*,err=9000)reduced_c_factor

  CALL write_line('FINISHED: read_reduced_c_factor',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error reading reduced c factor:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_reduced_c_factor
