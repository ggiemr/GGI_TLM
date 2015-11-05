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
! SUBROUTINE open_output_file_write
! SUBROUTINE open_output_file_read
! SUBROUTINE close_output_file
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
  
  if (write_info_file) then
    open(unit=info_file_unit,file=trim(problem_name)//info_file_extn)
  else
    open(unit=info_file_unit,file='/dev/null')
  end if
  
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
!     SUBROUTINE open_file
!
! DESCRIPTION
!     open_file:
!
!     open file with the provided unit, file extension and integer tag
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
!
! NAME
!     SUBROUTINE open_output_file_write
!
! DESCRIPTION
!     open_output_file_write:
!
!     open output file with the provided unit, file extension and compression flag
!     If compression is required then this process is achieved using a named pipe
!     linked to gzip
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/11/2015 CJS
!
!
SUBROUTINE open_output_file_write(file_unit,file,compression_flag)

IMPLICIT NONE

integer :: file_unit
character*(*)	:: file
logical	::compression_flag

! local variables

character(len=256)	:: filename
character(len=256)	:: gzfilename
character(len=256)	:: Instruction
integer			:: iostat

! START

  if (.NOT.compression_flag) then
  
    filename=trim(file)
    open(UNIT=file_unit,FILE=filename)
  
  else
! compression is required

    filename=trim(file)
    gzfilename=trim(filename)//'.gz'

! the name of the pipe is 'filename'
    call SYSTEM("rm -f "//trim(filename)//" ;mkfifo "//trim(filename))

! start a background gzip process with the pipe as input
! gzip writes to the filename with additional .gz extensiion

    Instruction="gzip -9 -c < "//trim(filename)//" > "//trim(gzfilename)//" &"
    call SYSTEM(Instruction)

! open a connection to the pipe for writing
    open(UNIT=file_unit,file=trim(filename),iostat=iostat,ACTION='WRITE')
    
  end if

  RETURN
  
END SUBROUTINE open_output_file_write
!
! NAME
!     SUBROUTINE open_output_file_read
!
! DESCRIPTION
!     open_output_file_read:
!
!     open output file with the provided unit, file extension and compression flag
!     If compression is required then this process is achieved using a named pipe
!     linked to gzip
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/11/2015 CJS
!
!
SUBROUTINE open_output_file_read(file_unit,file,compression_flag)

IMPLICIT NONE

integer :: file_unit
character*(*)	:: file
logical	::compression_flag

! local variables

character(len=256)	:: filename
character(len=256)	:: gzfilename
character(len=256)	:: Instruction
integer			:: iostat
logical			:: exists

! START

  if (.NOT.compression_flag) then
  
    filename=trim(file)
    open(UNIT=file_unit,FILE=filename)
  
  else
! compression is required

    filename=trim(file)
    gzfilename=trim(filename)//'.gz'

! check that the .gz file exists...
    INQUIRE(FILE=trim(gzfilename), EXIST=exists)
    if (.NOT.exists) then
      write(*,*)'File does not exist:',trim(gzfilename)
      STOP
    end if

! the name of the pipe is 'filename'
    call SYSTEM("rm -f "//trim(filename)//" ;mkfifo "//trim(filename))

! start a background gzip -dc process with the .gz file as input and the pipe as output

    Instruction="gzip -dc "//trim(gzfilename)//" > "//trim(filename)//"  &"
    call SYSTEM(Instruction)

! open a connection to the pipe for writing
    open(UNIT=file_unit,file=trim(filename),iostat=iostat,ACTION='READ')
    
  end if

  RETURN
  
END SUBROUTINE open_output_file_read
!
! NAME
!     SUBROUTINE close_output_file
!
! DESCRIPTION
!     close an output file used for reading or writing
!    
!
! COMMENTS
!     
!
! HISTORY
!
!     started 5/11/2015 CJS
!
!
SUBROUTINE close_output_file(file_unit,file,compression_flag)

IMPLICIT NONE

integer :: file_unit
character*(*)	:: file
logical	::compression_flag

! local variables

character(len=256)	:: filename
character(len=256)	:: gzfilename
character(len=256)	:: Instruction

! START

  if (.NOT.compression_flag) then

    close(unit=file_unit)
    
  else
  
    filename=trim(file)
    gzfilename=trim(filename)//'.gz'
  
! close the fortran connection to the pipe
    close(unit=file_unit)

! remove the named pipe
    Instruction="rm -f "//trim(filename)
    CALL system(Instruction)

  end if

  RETURN
  
END SUBROUTINE close_output_file
