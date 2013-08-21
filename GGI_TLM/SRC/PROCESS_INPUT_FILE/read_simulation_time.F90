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
! SUBROUTINE read_simulation_time
!
! NAME
!     read_simulation_time
!
! DESCRIPTION
!     read simulation time (in seconds)
!
! Example packet:
!
!Simulation_time
!1e-6   (real)
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE read_simulation_time

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: read_simulation_time',0,output_to_screen_flag)

  read(input_file_unit,*,err=9000)simulation_time

! check simulation time is greater than zero
  if ( simulation_time.lt.0d0 ) GOTO 9010

  CALL write_line('FINISHED: read_simulation_time',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error reading simulation_time packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading simulation_time packet',0,.TRUE.)
     CALL write_line('Simulation time should be greater than 0',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_simulation_time
