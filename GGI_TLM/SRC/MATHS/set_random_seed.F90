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
! SUBROUTINE set_random_seed
!
! NAME
!     set_random_seed
!
! DESCRIPTION
!     Set the random number generator seed value to a random number (always different)
!     Randomising the seed value is based on the subroutine on the webpage:
!
!     http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 29/01/2014 CJS
!
!
SUBROUTINE set_random_seed

IMPLICIT NONE

! local variables

  real*8  :: r_random
  
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt_local(8), pid, t(2), s
  integer(8) :: count, tms
  
! START  

  call random_seed(size = n)
    
  ALLOCATE(seed(n))
    
  
  ! XOR the current time and pid. The PID is
  ! useful in case one launches multiple instances of the same
  ! program in parallel.
  call system_clock(count)
  if (count /= 0) then
     t = transfer(count, t)
  else
     call date_and_time(values=dt_local)
     tms = (dt_local(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
  	  + dt_local(2) * 31_8 * 24 * 60 * 60 * 1000 &
  	  + dt_local(3) * 24 * 60 * 60 * 60 * 1000 &
  	  + dt_local(5) * 60 * 60 * 1000 &
  	  + dt_local(6) * 60 * 1000 + dt_local(7) * 1000 &
  	  + dt_local(8)
     t = transfer(tms, t)
  end if
  s = ieor(t(1), t(2))
  pid = getpid() + 1099279 ! Add a prime
  s = ieor(s, pid)
  if (n >= 3) then
     seed(1) = t(1) + 36269
     seed(2) = t(2) + 72551
     seed(3) = pid
     if (n > 3) then
  	seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
     end if
  else
     seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
  end if
 
  call random_seed(put=seed)

  DEALLOCATE(seed)

! evaluate a few random numbers to 'warm up' the random number generator if small values of seed are used for example  
  do i=1,256
    CALL random_number(r_random)
  end do
       
  RETURN
       
END SUBROUTINE set_random_seed
