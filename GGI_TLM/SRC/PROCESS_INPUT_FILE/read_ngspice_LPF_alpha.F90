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
! SUBROUTINE read_ngspice_LPF_alpha
!
! NAME
!     read_ngspice_LPF_alpha
!
! DESCRIPTION
!     read the ngspice low pass filter coefficient, alpha
!
! Example packet:
!
!
!ngspice_LPF_alpha
!0.25
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/04/2019 CJS Introduce ngspice link
!
!
SUBROUTINE read_ngspice_LPF_alpha

USE TLM_general
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: read_ngspice_LPF_alpha',0,output_to_screen_flag)

  read(input_file_unit,*,err=9000)LPF_alpha
  
  if (LPF_alpha.LT.0.0) then
    CALL write_line('low pass filter coefficient should be greater than zero',0,output_to_screen_flag)
    GOTO 9000
  end if

  CALL write_line('FINISHED: read_ngspice_LPF_alpha',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error reading low pass filter coefficient:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
    
END SUBROUTINE read_ngspice_LPF_alpha
