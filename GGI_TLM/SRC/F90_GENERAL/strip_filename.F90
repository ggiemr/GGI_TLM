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
!
! NAME
!     strip_path
!
! DESCRIPTION
!     takes a filename with a full path and returns only the filename
!
! HISTORY
!
!     started 7/10/12 CJS
!
! COMMENTS
!     
FUNCTION strip_path(input_string) RESULT(result_string)

USE File_general

IMPLICIT NONE

! variables passed to subroutine

  character(LEN=filename_length) :: input_string
  
! result  
  character(LEN=filename_length) :: result_string
  
! local variables

  integer i
  integer last_slash
  integer full_length
  integer length
  
! START
  
  full_length=len_trim(input_string)

  last_slash=0
  result_string=''
  
  do i=1,full_length
  
    if (input_string(i:i).eq.'/') last_slash=i
    
  end do
  
  length=full_length-last_slash
    
  result_string(1:length)=input_string(last_slash+1:last_slash+length)
  
  return
  
end FUNCTION strip_path
