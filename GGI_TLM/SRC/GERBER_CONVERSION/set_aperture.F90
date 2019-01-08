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
!     SUBROUTINE set_aperture
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
SUBROUTINE set_aperture(line,pos,nDcodes,maxDcodes,Dcode,aperture)

IMPLICIT NONE

character*256 :: line
integer :: pos

integer :: nDcodes
integer :: maxDcodes
integer :: Dcode(1:maxDcodes)
integer :: aperture

integer :: dc
integer :: i

character :: eol

! START

CALL read_integer_from_line(line,pos,dc)

! check we have a valid Dcode
if (dc.LT.10) then
  write(*,*)'Invalid Dcode. Dcode is less than 10'
  STOP
end if

! look through the Dcode list for the specified Dcode
aperture=0

do i=1,nDcodes

  if (Dcode(i).EQ.dc) then
    aperture=i
    EXIT
  end if
  
end do ! next dcode

if (aperture.EQ.0) then
  write(*,*)'Dcode not found: Dcode=',dc
end if

write(*,*)'Setting aperture D',dc,' position in aperture list=',aperture

RETURN

END SUBROUTINE set_aperture
