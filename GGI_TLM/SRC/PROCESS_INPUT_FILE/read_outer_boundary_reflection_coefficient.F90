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
! SUBROUTINE read_outer_boundary_reflection_coefficient
!
! NAME
!     read_outer_boundary_reflection_coefficient
!
! DESCRIPTION
!     read outer boundary reflection coefficients
!
! Example packet:
!
!Outer_boundary_reflection_coefficient
! R_xmin R_xmax   R_ymin R_ymax   R_zmin R_zmax    (6*real)
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE read_outer_boundary_reflection_coefficient

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: Read_outer_boundary_reflection_coefficients',0,output_to_screen_flag)

! read reflection coefficient data
  read(input_file_unit,*,err=9000)R_xmin,R_xmax,R_ymin,R_ymax,R_zmin,R_zmax

! check to ensure that reflection coefficients are less than or equal to 1
  if ( (abs(R_xmin).GT.1d0).OR.(abs(R_xmax).GT.1d0).OR.	&
       (abs(R_ymin).GT.1d0).OR.(abs(R_ymax).GT.1d0).OR.	&
       (abs(R_zmin).GT.1d0).OR.(abs(R_zmax).GT.1d0) ) GOTO 9010

  CALL write_line('FINISHED: Read_outer_boundary_reflection_coefficients',0,output_to_screen_flag)
       
  RETURN
  
9000 CALL write_line('Error reading outer_boundary_reflection_coefficient packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error magnitude of outer boundary reflection coefficients should be lsee than or equal to zero',0,.TRUE.)
     CALL write_line('input_filename',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_outer_boundary_reflection_coefficient
