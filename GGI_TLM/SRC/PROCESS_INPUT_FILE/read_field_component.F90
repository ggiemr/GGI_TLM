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
! SUBROUTINE read_field_component
! SUBROUTINE read_field_component_surface_output
!
! NAME
!     read_field_component
!
! DESCRIPTION
!     read a field component
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE read_field_component(file_unit,field_component)

USE TLM_general

IMPLICIT NONE

  integer		:: file_unit
  integer		:: field_component

! local variables

character*256	:: input_line

! START  
    read(file_unit,'(A)')input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,256)

    if (input_line.EQ.Ex_string) then
      field_component=Ex
    else if (input_line.EQ.Ey_string) then
      field_component=Ey
    else if (input_line.EQ.Ez_string) then
      field_component=Ez
    else if (input_line.EQ.Hx_string) then
      field_component=Hx
    else if (input_line.EQ.Hy_string) then
      field_component=Hy
    else if (input_line.EQ.Hz_string) then
      field_component=Hz
    else
      GOTO 9010
    end if  
    
  RETURN
  
9000 CALL write_line('Error reading field component from file:',0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
     
9010 CALL write_line('Error reading field component',0,.TRUE.)
     CALL write_line("Expecting field compoent 'Ex', 'Ey', 'Ez', 'Hx', 'Hy' or 'Hz'",0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
  
END SUBROUTINE read_field_component
! SUBROUTINE read_field_component_surface_output
!
! NAME
!     read_field_component_surface_output
!
! DESCRIPTION
!     read a field component for surface output
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 17/01/2013 CJS
!
!
SUBROUTINE read_field_component_surface_output(file_unit,field_component)

USE TLM_general

IMPLICIT NONE

  integer		:: file_unit
  integer		:: field_component

! local variables

character*2	:: input_line

! START  
    read(file_unit,'(A2)')input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,2)

    if (input_line.EQ.Ex_string) then
      field_component=Ex
    else if (input_line.EQ.Ey_string) then
      field_component=Ey
    else if (input_line.EQ.Ez_string) then
      field_component=Ez
    else if (input_line.EQ.Hx_string) then
      field_component=Hx
    else if (input_line.EQ.Hy_string) then
      field_component=Hy
    else if (input_line.EQ.Hz_string) then
      field_component=Hz
    else if (input_line.EQ.Jx_string) then
      field_component=Jx
    else if (input_line.EQ.Jy_string) then
      field_component=Jy
    else if (input_line.EQ.Jz_string) then
      field_component=Jz
    else if (input_line.EQ.Emagnitude_string) then
      field_component=Emagnitude
    else if (input_line.EQ.Hmagnitude_string) then
      field_component=Hmagnitude
    else if (input_line.EQ.Jmagnitude_string) then
      field_component=Jmagnitude
    else
      GOTO 9010
    end if  
    
  RETURN
  
9000 CALL write_line('Error reading field component for surface_output from file:',0,.TRUE.)
     CALL write_line('input_filename',0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
     
9010 CALL write_line('Error reading field component for surface_output',0,.TRUE.)
     CALL write_line("Expecting field compoent 'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', 'Em', 'Hm', 'Jx', 'Jy', 'Jz', 'Jm', 'Power'"&
                     ,0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
  
END SUBROUTINE read_field_component_surface_output
! SUBROUTINE read_field_component_frequency_surface_output
!
! NAME
!     read_field_component_frequency_surface_output
!
! DESCRIPTION
!     read a field component for surface output
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 17/01/2013 CJS
!
!
SUBROUTINE read_field_component_frequency_surface_output(file_unit,field_component)

USE TLM_general

IMPLICIT NONE

  integer		:: file_unit
  integer		:: field_component

! local variables

character*2	:: input_line

! START  
    read(file_unit,'(A2)')input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,2)

    if (input_line.EQ.Ex_string) then
      field_component=Ex
    else if (input_line.EQ.Ey_string) then
      field_component=Ey
    else if (input_line.EQ.Ez_string) then
      field_component=Ez
    else if (input_line.EQ.Hx_string) then
      field_component=Hx
    else if (input_line.EQ.Hy_string) then
      field_component=Hy
    else if (input_line.EQ.Hz_string) then
      field_component=Hz
    else if (input_line.EQ.Jx_string) then
      field_component=Jx
    else if (input_line.EQ.Jy_string) then
      field_component=Jy
    else if (input_line.EQ.Jz_string) then
      field_component=Jz
    else
      GOTO 9010
    end if  
    
  RETURN
  
9000 CALL write_line('Error reading field component for surface_output from file:',0,.TRUE.)
     CALL write_line('input_filename',0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
     
9010 CALL write_line('Error reading field component for surface_output',0,.TRUE.)
     CALL write_line("Expecting field compoent 'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', 'Em', 'Hm', 'Jx', 'Jy', 'Jz', 'Jm', 'Power'"&
                     ,0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
  
END SUBROUTINE read_field_component_frequency_surface_output
