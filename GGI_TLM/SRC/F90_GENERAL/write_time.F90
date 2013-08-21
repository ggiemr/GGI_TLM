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
!	write_date_and_time
!
! DESCRIPTION
!     write the date and time to the given file unit
!
! HISTORY
!
!     started 10/02/09 CJS
!
! COMMENTS
!     
subroutine write_date_and_time(unit)

! input variables
integer unit

! local variables

character*8 date
character*10 time
character*5 zone
integer values(8)

integer i

  call DATE_AND_TIME(date,time,zone,values)
  
  if (unit.ne.0) then
    write(unit,*)'Date:',date(7:8),'/',date(5:6),'/',date(1:4)
    write(unit,*)'Time:',time(1:2),':',time(3:4),':',time(5:10)
  else
    write(*,*)'Date:',date(7:8),'/',date(5:6),'/',date(1:4)
    write(*,*)'Time:',time(1:2),':',time(3:4),':',time(5:10)
  end if

return
end
