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
!SUBROUTINE write_license
!SUBROUTINE write_progress
!SUBROUTINE write_line
!SUBROUTINE write_line_integer
!SUBROUTINE write_line_real
!
! NAME
!	write_license
!
! DESCRIPTION
!     writes the license agreement note 
!
! HISTORY
!
!     started 3/05/13 CJS
!
! COMMENTS
!     
SUBROUTINE write_license()

! variables passed to subroutine
  
! local variables
  
! START

  write(*,*)'  '
  write(*,*)'  Copyright (C) 2013  Chris Smartt'
  write(*,*)'  This program comes with ABSOLUTELY NO WARRANTY'
  write(*,*)'  This is free software, and you are welcome to redistribute it '
  write(*,*)'  under certain conditions.'
  write(*,*)'  See the GNU General Public License for more details.'
  write(*,*)'  You should have received a copy of the GNU General Public License'
  write(*,*)'  along with this program.  If not, see <http://www.gnu.org/licenses/>.'
  write(*,*)'  '

  
  RETURN
  
END SUBROUTINE write_license
!
! NAME
!	write_progress
!
! DESCRIPTION
!     writes the given string to the progress file
!
! HISTORY
!
!     started 12/10/12 CJS
!
! COMMENTS
!     
SUBROUTINE write_progress(line)

  USE File_information

! variables passed to subroutine

  character*(*)	:: line
  
! local variables
  
! START
  
  open(unit=progress_file_unit,file=progress_filename)
  
  write(progress_file_unit,'(A)')trim(line)
  
  close(unit=progress_file_unit)
  
  RETURN
  
END SUBROUTINE write_progress
!
! NAME
!	write_line
!
! DESCRIPTION
!     writes the given string to the desired file unit or screen 
!
! HISTORY
!
!     started 7/08/12 CJS
!
! COMMENTS
!     
SUBROUTINE write_line(line,unit,flag)

! variables passed to subroutine

  character*(*)	:: line
  integer	:: unit
  logical	:: flag
  
! local variables
  
! START

  if (.not.flag) RETURN
  
  if (unit.eq.0) then
    write(*,'(A)')trim(line)
    flush(5)
  else
    write(unit,'(A)')trim(line)
    flush(unit)
  end if
  
  RETURN
  
END SUBROUTINE write_line
!
! NAME
!	write_line_integer
!
! DESCRIPTION
!     writes the given string and integer to the desired file unit or screen 
!
! HISTORY
!
!     started 7/08/12 CJS
!
! COMMENTS
!     
SUBROUTINE write_line_integer(line,op_int,unit,flag)

! variables passed to subroutine

  character*(*)	:: line
  integer	:: op_int
  integer	:: unit
  logical	:: flag
  
! local variables
  
! START

  if (.not.flag) RETURN
  
  if (unit.eq.0) then
    write(*,'(A,I10)')trim(line),op_int
    flush(5)
  else
    write(unit,'(A,I10)')trim(line),op_int
    flush(unit)
  end if
  
  RETURN
  
END SUBROUTINE write_line_integer
!
! NAME
!	write_line_real
!
! DESCRIPTION
!     writes the given string and real number to the desired file unit or screen 
!
! HISTORY
!
!     started 7/08/12 CJS
!
! COMMENTS
!     
SUBROUTINE write_line_real(line,op_real,unit,flag)

! variables passed to subroutine

  character*(*)	:: line
  real*8	:: op_real
  integer	:: unit
  logical	:: flag
  
! local variables
  
! START

  if (.not.flag) RETURN
  
  if (unit.eq.0) then
    write(*,'(A,E16.6)')trim(line),op_real
    flush(5)
  else
    write(unit,'(A,E16.6)')trim(line),op_real
    flush(unit)
  end if
  
  
  RETURN
  
END SUBROUTINE write_line_real
