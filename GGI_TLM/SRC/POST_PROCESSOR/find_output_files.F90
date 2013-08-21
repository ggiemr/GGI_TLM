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
! SUBROUTINE find_output_files
!
! NAME
!    find_output_files
!
! DESCRIPTION
!     
!     Find all the output files generated with the given problem name
!
!  Output data supported:
!     1. Time Domain excitation function
!     2. Time Domain field output
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
SUBROUTINE find_output_files

USE TLM_general
USE file_information

IMPLICIT NONE

! local variables

character(len=256)	:: filename
logical			:: file_exists

! START

  
  CALL write_line('CALLED: find_output_files',0,output_to_screen_flag)

! Check for Time Domain field excitation function file  
  
  filename=trim(problem_name)//excitation_output_extn 
  inquire(file=trim(filename),exist=file_exists)
  if( file_exists ) then
    write(*,*)
    write(*,*)'Time domain excitation function output file:'
    write(*,*)trim(filename)
  end if

! Check for Time Domain field output file  

  filename=trim(problem_name)//field_output_extn
  inquire(file=trim(filename),exist=file_exists)
  if( file_exists ) then
    write(*,*)
    write(*,*)'Time domain field output file:'
    write(*,*)trim(filename)
  end if
  
  write(*,*)
  CALL write_line('FINISHED: find_output_files',0,output_to_screen_flag)

  RETURN
  
  
END SUBROUTINE find_output_files
