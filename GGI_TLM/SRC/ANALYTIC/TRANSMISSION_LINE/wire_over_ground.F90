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
! Name wire_over_ground
!     
!
! Description
!     calculate the inductance, capacitance and characteristic impedance of a wire over a ground plane
!
! Comments:
!      
!
! History
!
!     started 7/3/2018 CJS
!

PROGRAM wire_over_ground

USE constants

IMPLICIT NONE

! local_variables

  real r,h
  real t,s
  real Z_wire_over_ground
  real C_wire_over_ground
  real L_wire_over_ground
   
! function_types

! START

!  write(*,*)'CALLED: wire_over_ground'

  write(*,*)'Enter wire radius (mm)'
  read(*,*)r
  write(*,*)'Enter height over ground plane (mm)'
  read(*,*)h
  
  s=2.0*h
  t=s/(2.0*r)
  
  C_wire_over_ground=2.0*pi*eps0/(log(t+sqrt(t*t-1.0)))   ! note this is twice C for the two wire case
  
  L_wire_over_ground=mu0*eps0/C_wire_over_ground
  
  Z_wire_over_ground=sqrt(L_wire_over_ground/C_wire_over_ground)
    
  write(*,*)'Capacitance per unit length=',C_wire_over_ground,' F/m'
  
  write(*,*)'Inductance per unit length =',L_wire_over_ground,' H/m'
   
  write(*,*)'Characteristic Impedance=',Z_wire_over_ground,' ohms'
   
  write(*,*)'Velocity=',1.0/sqrt(L_wire_over_ground*C_wire_over_ground),' m/s'
 
!  write(*,*)'FINISHED: wire_over_ground'
  
END 

  
  
  
  
