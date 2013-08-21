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
  SUBROUTINE calc_phi_const(rp,tp,rw,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  real*8 phi_local

! START

  if (rp.lt.rw) then
    call calc_phi_ICR_const(rp,tp,rw,phi_local)
  else
    call calc_phi_OCR_const(rp,tp,rw,phi_local)
  end if
 
 return
 
 END SUBROUTINE calc_phi_const
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_cos(rp,tp,rw,m,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  integer m
  real*8 phi_local

! START
 
  if (rp.lt.rw) then
    call calc_phi_ICR_cos(rp,tp,rw,m,phi_local)
  else
    call calc_phi_OCR_cos(rp,tp,rw,m,phi_local)
  end if
 
 return
 
 END SUBROUTINE calc_phi_cos
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_sin(rp,tp,rw,m,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  integer m
  real*8 phi_local

! START
 
  if (rp.lt.rw) then
    call calc_phi_ICR_sin(rp,tp,rw,m,phi_local)
  else
    call calc_phi_OCR_sin(rp,tp,rw,m,phi_local)
  end if
 
 return
 
 END SUBROUTINE calc_phi_sin
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_ICR_const(rp,tp,rw,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  real*8 phi_local

! START
 
 phi_local=-rw*log(rw)/eps0
 
 return
 
 END SUBROUTINE calc_phi_ICR_const
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_ICR_cos(rp,tp,rw,m,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  integer m
  real*8 phi_local

! START
 
 phi_local=(rp**m)*cos(m*tp)/( 2d0*eps0*m*(rw**(m-1)) )
 
 return
 
 END SUBROUTINE calc_phi_ICR_cos
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_ICR_sin(rp,tp,rw,m,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  integer m
  real*8 phi_local

! START
 
 phi_local=(rp**m)*sin(m*tp)/( 2d0*eps0*m*(rw**(m-1)) ) 
 
 return
 
 END SUBROUTINE calc_phi_ICR_sin
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_OCR_const(rp,tp,rw,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  real*8 phi_local

! START
 
 phi_local=-rw*log(rp)/eps0
 
 return
 
 END SUBROUTINE calc_phi_OCR_const
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_OCR_cos(rp,tp,rw,m,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  integer m
  real*8 phi_local

! START
 
 phi_local=(rw**(m+1))*cos(m*tp)/( 2d0*eps0*m*(rp**m) ) 
 
 return
 
 END SUBROUTINE calc_phi_OCR_cos
!
! _________________________________________________
!
!
  SUBROUTINE calc_phi_OCR_sin(rp,tp,rw,m,phi_local)
   
USE constants

IMPLICIT NONE
  
! variables passed to subroutine 
  
  real*8 rw,rp,tp
  integer m
  real*8 phi_local

! START
 
 phi_local=(rw**(m+1))*sin(m*tp)/( 2d0*eps0*m*(rp**m) )
 
 return
 
 END SUBROUTINE calc_phi_OCR_sin
!
! _________________________________________________
!
!

