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
! SUBROUTINE TLM_set_small_to_zero
!
! NAME
!     TLM_set_small_to_zero
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/11/2020 CJS try to eliminate a problem which may be due to denormal numbers
!
!
SUBROUTINE TLM_set_small_to_zero

USE TLM_general
USE mesh

IMPLICIT NONE

! local variables

  integer cx,cy,cz
  integer i

! START
  
  CALL write_line('CALLED: TLM_set_small_to_zero',0,timestepping_output_to_screen_flag)
    
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
      
        do i=1,12
	
          if (abs(V(i,cx,cy,cz)).LT.small_to_zero_value) V(i,cx,cy,cz)=0d0
	
        end do  ! next port
	
      end do  ! next z cell
    end do    ! next y cell
  end do      ! next x cell
  
  CALL write_line('FINISHED: TLM_set_small_to_zero',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE TLM_set_small_to_zero
