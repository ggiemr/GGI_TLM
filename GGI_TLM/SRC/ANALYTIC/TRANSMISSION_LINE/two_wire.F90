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
! Name two_wire
!     
!
! Description
!     calculate the inductance, capacitance and characteristic impedance of a two wire transmission line
!
! Comments:
!      
!
! History
!
!     started 7/3/2018 CJS
!

PROGRAM two_wire

USE constants

IMPLICIT NONE

! local_variables

  real r,s
  real t
  real Z_two_wire
  real C_two_wire
  real L_two_wire
   
! function_types

! START

!  write(*,*)'CALLED: two_wire'

  write(*,*)'Enter wire radius (mm)'
  read(*,*)r
  write(*,*)'Enter the wire separation (mm)'
  read(*,*)s
  
  t=s/(2.0*r)
  
  C_two_wire=pi*eps0/(log(t+sqrt(t*t-1.0)))  
  
  L_two_wire=mu0*eps0/C_two_wire
  
  Z_two_wire=sqrt(L_two_wire/C_two_wire)
    
  write(*,*)'Capacitance per unit length=',C_two_wire,' F/m'
  
  write(*,*)'Inductance per unit length =',L_two_wire,' H/m'
   
  write(*,*)'Characteristic Impedance=',Z_two_wire,' ohms'
   
  write(*,*)'Velocity=',1.0/sqrt(L_two_wire*C_two_wire),' m/s'
 
!  write(*,*)'FINISHED: two_wire'
  
END 

  
  
  
  
