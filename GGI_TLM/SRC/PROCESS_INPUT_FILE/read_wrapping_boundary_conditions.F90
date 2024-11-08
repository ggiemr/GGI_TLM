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
! SUBROUTINE read_wrapping_boundary_conditions
!
! NAME
!     read_wrapping_boundary_conditions
!
! DESCRIPTION
!     read wrapping boundary conditions
!
! Example packet:
!
!
!Wrapping_boundary_conditions
!1 1 0
!
! COMMENTS
!     
!
! HISTORY
!
!     started 25/04/2013 CJS
!
!
SUBROUTINE read_wrapping_boundary_conditions

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

  integer	:: WR_x,WR_y,WR_z
  
! START  

  CALL write_line('CALLED: Read_wrapping_boundary_conditions',0,output_to_screen_flag)

! read reflection coefficient data
  read(input_file_unit,*,err=9000,end=9000)WR_x,WR_y,WR_z

! check to ensure that wrapping boundary coefficients are equal to 0 or 1
  if ((WR_x.NE.1).AND.(WR_x.NE.0)) GOTO 9010
  if ((WR_y.NE.1).AND.(WR_y.NE.0)) GOTO 9010
  if ((WR_z.NE.1).AND.(WR_z.NE.0)) GOTO 9010
  
  if (WR_x.EQ.1) wrap_x=.TRUE.
  if (WR_y.EQ.1) wrap_y=.TRUE.
  if (WR_z.EQ.1) wrap_z=.TRUE.

  CALL write_line('FINISHED: Read_wrapping_boundary_conditionss',0,output_to_screen_flag)
       
  RETURN
  
9000 CALL write_line('Error reading wrapping_boundary_conditions packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
  
9010 CALL write_line('Error reading wrapping_boundary_conditions packet data from input file:',0,.TRUE.)
     CALL write_line('Wrapping boundary falg should be 1 or 0 ',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
       
END SUBROUTINE read_wrapping_boundary_conditions
