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
  SUBROUTINE write_matrix(M,n,scale,unit)

IMPLICIT NONE

! variables passed to subroutine

  integer n
  real*8 :: M(1:n,1:n)
  real*8 scale
  integer unit
  
! local variables 
   
  integer row,col

! START

  if (n.le.10) then  ! write as a matrix
 
    do row=1,n
  
        if (unit.eq.0) then
	  write(*,8000)(M(row,col)*scale,col=1,n)
	else
	  write(unit,8000)(M(row,col)*scale,col=1,n)
	end if  
8000    format(10F9.4)
 
    end do ! next row
    
  else  ! too large so write in rows
  
    do row=1,n
  
      do col=1,n
 
        if (unit.eq.0) then
	  write(*,8010)row,col,M(row,col)*scale
	else
	  write(unit,8010)row,col,M(row,col)*scale
	end if  
8010    format(2I6,E16.5)
 
      end do ! next col
    
    end do ! next row
  
  end if
  
  return
  
  end SUBROUTINE write_matrix
!  
! ________________________________________________  
!  
!  
  SUBROUTINE write_vector(M,n,scale,unit)

IMPLICIT NONE

! variables passed to subroutine

  integer n
  real*8 :: M(1:n)
  real*8 scale
  integer unit
  
! local variables 
   
  integer row,col

! START

  
  do row=1,n
  
    if (unit.eq.0) then
      write(*,8010)row,M(row)*scale
    else
      write(unit,8010)row,M(row)*scale
    end if  
8010    format(I6,E16.5)
 
  end do ! next row
  
  return
  
  end SUBROUTINE write_vector
!  
! ________________________________________________  
!  
!  
  SUBROUTINE write_matrix_noscale(M,n,unit)

IMPLICIT NONE

! variables passed to subroutine

  integer n
  real*8 :: M(1:n,1:n)
  integer unit
  
! local variables 
   
  integer row,col

! START

  
  do row=1,n
  
    do col=1,n
 
      if (unit.eq.0) then
        write(*,8010)row,col,M(row,col)
      else
	write(unit,8010)row,col,M(row,col)
      end if  
8010  format(2I6,E16.5)
 
    end do ! next col
    
  end do ! next row
  
  
  return
  
  end SUBROUTINE write_matrix_noscale
