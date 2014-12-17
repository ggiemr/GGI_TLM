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
PROGRAM GGI_TLM_generate_random_numbers 

real*8  in1,in2

real*8  result

integer n

integer i

character ch
character*256 filename

! START

write(*,*)'Enter the number of random real numbers to generate'

read(*,*)n

write(*,*)'Enter the minimum and maximum possible values for the random number'

read(*,*)in1,in2

write(*,*)'Do you want to randomise the random number generator seed values (y/n)'
read(*,'(A)')ch

if ((ch.eq.'y').OR.(ch.eq.'Y')) then
  CALL set_random_seed()
end if

write(*,*)'Enter the filename for the random number data'

read(*,'(A)')filename

open(unit=10,file=trim(filename))

do i=1,n

  CALL random_number(result)
  write(*,*)in1+(in2-in1)*result
  write(10,*)in1+(in2-in1)*result

end do

close(unit=10)

END PROGRAM GGI_TLM_generate_random_numbers 
