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
! SUBROUTINE open_general_files
! SUBROUTINE close_general_files
! SUBROUTINE open_vtk_file
! SUBROUTINE close_vtk_file
! SUBROUTINE open_file
! SUBROUTINE close_file
! NAME
!     SUBROUTINE open_general_files
!
! DESCRIPTION
!     open_general_files:
!
!     read the problem name and open the following files:
!     input file
!     warning file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE open_general_files

USE TLM_general
USE file_information

IMPLICIT NONE

! local variables

! START

  
  CALL write_line('CALLED: open_general_files',0,output_to_screen_flag)
  
  open(unit=input_file_unit,file=trim(problem_name)//input_file_extension,status='OLD',err=9000)
  
  open(unit=info_file_unit,file=trim(problem_name)//info_file_extn)
  
  open(unit=warning_file_unit,file=trim(problem_name)//warning_file_extension)
  
  CALL write_line('FINISHED: open_general_files',0,output_to_screen_flag)

  RETURN
  
9000 CALL write_line('Error opening input file:',0,.TRUE.)
     CALL write_line(trim(problem_name)//input_file_extension,0,.TRUE.)
     STOP
  
END SUBROUTINE open_general_files
!
! NAME
!     SUBROUTINE close_general_files
!
! DESCRIPTION
!     close_general_files:
!
!     Close the following files:
!     input file
!     warning file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE close_general_files

USE file_information
USE TLM_general

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: close_general_files',0,output_to_screen_flag)
  
  close(unit=input_file_unit)
  
  close(unit=info_file_unit)
  
  close(unit=warning_file_unit)
  
  CALL write_line('FINISHED: close_general_files',0,output_to_screen_flag)

  RETURN
  
END SUBROUTINE close_general_files
!
! NAME
!     SUBROUTINE open_vtk_file
!
! DESCRIPTION
!     open_vtk_file:
!
!     open vtk file with the provided unit, file extension and integer tag
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE open_vtk_file(file_unit,file_extension,number)

USE TLM_general
USE file_information

IMPLICIT NONE

integer :: file_unit
character*(*)	:: file_extension
integer	:: number

! local variables

character(len=256)	:: filename
character(len=256)	:: temp_filename

! START

    temp_filename=trim(problem_name)//trim(file_extension)
    CALL add_integer_to_filename(temp_filename,number,filename)
    
    open(UNIT=file_unit,FILE=filename)

  RETURN
  
END SUBROUTINE open_vtk_file
!
! NAME
!     SUBROUTINE close_vtk_file
!
! DESCRIPTION
!     close_vtk_file:
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE close_vtk_file(file_unit)

USE TLM_general
USE file_information

IMPLICIT NONE

integer :: file_unit

! local variables

! START

  close(unit=file_unit)

  RETURN
  
END SUBROUTINE close_vtk_file
!
! NAME
!     SUBROUTINE open_vtk_file
!
! DESCRIPTION
!     open_vtk_file:
!
!     open vtk file with the provided unit, file extension and integer tag
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE open_file(file_unit,file_extension)

USE TLM_general
USE file_information

IMPLICIT NONE

integer :: file_unit
character*(*)	:: file_extension

! local variables

character(len=256)	:: filename

! START

    filename=trim(problem_name)//trim(file_extension)
    
    open(UNIT=file_unit,FILE=filename)

  RETURN
  
END SUBROUTINE open_file
!
! NAME
!     SUBROUTINE close_file
!
! DESCRIPTION
!     close_file:
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE close_file(file_unit)

USE TLM_general
USE file_information

IMPLICIT NONE

integer :: file_unit

! local variables

! START

  close(unit=file_unit)

  RETURN
  
END SUBROUTINE close_file
