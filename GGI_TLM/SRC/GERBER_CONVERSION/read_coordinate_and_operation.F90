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
!     SUBROUTINE read_coordinate_and_operation
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 10/12/19 CJS
!     
!
SUBROUTINE read_coordinate_and_operation(line,x,y,xfni,xfnd,yfni,yfnd,s,operation)
        
IMPLICIT NONE

character*256 :: line
real*8  :: x,y
integer :: xfni,xfnd,yfni,yfnd
real*8  :: s
integer :: operation

integer :: pos

logical :: have_coord
logical :: have_operation

! START

pos=1
have_coord=.FALSE.
have_operation=.FALSE.

DO

  if (line(pos:pos).EQ.'X') then

    pos=pos+1
    CALL read_formatted_coordinate_from_line(line,pos,x,xfni,xfnd,s)
    have_coord=.TRUE.
  
  else if (line(pos:pos).EQ.'Y') then

    pos=pos+1
    CALL read_formatted_coordinate_from_line(line,pos,y,yfni,yfnd,s)
    have_coord=.TRUE.
    
  else if (line(pos:pos).EQ.'D') then

    pos=pos+1
    CALL read_integer_from_line(line,pos,operation)
    have_operation=.TRUE.
  
  else if (line(pos:pos).EQ.'*') then

    write(*,*)'New coordinate:',x,y,' operation=',operation
    if (have_coord.AND.have_operation) then
      RETURN
    else
      write(*,*)'We are missing either a coordinate or operation number'
      STOP
    end if
  
  else

    write(*,*)'Error reading coordinates: found character:',line(pos:pos)
    STOP
  
  end if

END DO

END SUBROUTINE read_coordinate_and_operation
