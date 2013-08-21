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
! SUBROUTINE write_error_line
!
! NAME
!     write_error_line
!
! DESCRIPTION
!     
!     write the line of an input file which triggered an error in
!     reading the input file
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE write_error_line(file_unit)

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

integer	:: file_unit

! local variables

character*256	:: input_line

! START  

  CALL write_line('Error line:',0,.TRUE.)
  
  backspace(unit=file_unit)
  
! read error line from input file
  read(file_unit,'(A)')input_line

  CALL write_line(trim(input_line),0,.TRUE.)
  
  RETURN
  
END SUBROUTINE write_error_line
