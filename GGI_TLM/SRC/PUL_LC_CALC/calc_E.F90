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
  SUBROUTINE calc_E_const(rp,tp,rb,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  real*8 Er_local,Et_local

! START

  if (rp.lt.rb) then
    call calc_E_ICR_const(rp,tp,rb,Er_local,Et_local)
  else
    call calc_E_OCR_const(rp,tp,rb,Er_local,Et_local)
  end if
 
 return
 
 END SUBROUTINE calc_E_const
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_cos(rp,tp,rb,m,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  integer m
  real*8 Er_local,Et_local

! START
 
  if (rp.lt.rb) then
    call calc_E_ICR_cos(rp,tp,rb,m,Er_local,Et_local)
  else
    call calc_E_OCR_cos(rp,tp,rb,m,Er_local,Et_local)
  end if
 
 return
 
 END SUBROUTINE calc_E_cos
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_sin(rp,tp,rb,m,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  integer m
  real*8 Er_local,Et_local

! START
 
  if (rp.lt.rb) then
    call calc_E_ICR_sin(rp,tp,rb,m,Er_local,Et_local)
  else
    call calc_E_OCR_sin(rp,tp,rb,m,Er_local,Et_local)
  end if
 
 return
 
 END SUBROUTINE calc_E_sin
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_ICR_const(rp,tp,rb,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  real*8 Er_local,Et_local

! START
 
 Er_local=0d0
 Et_local=0d0
 
 return
 
 END SUBROUTINE calc_E_ICR_const
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_ICR_cos(rp,tp,rb,m,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  integer m
  real*8 Er_local,Et_local

! START

  Er_local=-((rp/rb)**(m-1))*cos(m*tp)/( 2d0*eps0 ) 
  Et_local= ((rp/rb)**(m-1))*sin(m*tp)/( 2d0*eps0 )
 
 return
 
 END SUBROUTINE calc_E_ICR_cos
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_ICR_sin(rp,tp,rb,m,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  integer m
  real*8 Er_local,Et_local

! START
 
  Er_local=-((rp/rb)**(m-1))*sin(m*tp)/( 2d0*eps0 )
  Et_local=-((rp/rb)**(m-1))*cos(m*tp)/( 2d0*eps0 )
 
 return
 
 END SUBROUTINE calc_E_ICR_sin
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_OCR_const(rp,tp,rb,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  real*8 Er_local,Et_local

! START
 
  Er_local=rb/(eps0*rp)
  Et_local=0d0
 
 return
 
 END SUBROUTINE calc_E_OCR_const
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_OCR_cos(rp,tp,rb,m,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  integer m
  real*8 Er_local,Et_local

! START
 
  Er_local= ((rb/rp)**(m+1))*cos(m*tp)/( 2d0*eps0 )
  Et_local= ((rb/rp)**(m+1))*sin(m*tp)/( 2d0*eps0 )
 
 return
 
 END SUBROUTINE calc_E_OCR_cos
!
! _________________________________________________
!
!
  SUBROUTINE calc_E_OCR_sin(rp,tp,rb,m,Er_local,Et_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rb,rp,tp
  integer m
  real*8 Er_local,Et_local

! START
 
  Er_local= ((rb/rp)**(m+1))*sin(m*tp)/( 2d0*eps0 )
  Et_local=-((rb/rp)**(m+1))*cos(m*tp)/( 2d0*eps0 )
 
 return
 
 END SUBROUTINE calc_E_OCR_sin
!
! _________________________________________________
!
!

