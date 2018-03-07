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
! Name coax
!     
!
! Description
!     calculate the inductance, capacitance and characteristic impedance of a coaxial transmission line
!
! Comments:
!      
!
! History
!
!     started 7/3/2018 CJS
!

PROGRAM coax

USE constants

IMPLICIT NONE

! local_variables

  real rw,rs,epsr
  real Z_coax
  real C_coax
  real L_coax
   
! function_types

! START

!  write(*,*)'CALLED: coax'

  write(*,*)'Enter inner wire radius (mm)'
  read(*,*)rw
  write(*,*)'Enter shield radius (mm)'
  read(*,*)rs
  write(*,*)'Enter the dielectric relative permittivity'
  read(*,*)epsr
  
  C_coax=2.0*pi*eps0*epsr/log(rs/rw)
  
  L_coax=(mu0/(2.0*pi))*log(rs/rw)
  
  Z_coax=sqrt(L_coax/C_coax)
    
  write(*,*)'Capacitance per unit length=',C_coax,' F/m'
  
  write(*,*)'Inductance per unit length =',L_coax,' H/m'
   
  write(*,*)'Characteristic Impedance=',Z_coax,' ohms'
   
  write(*,*)'Velocity=',1.0/sqrt(L_coax*C_coax),' m/s'
 
!  write(*,*)'FINISHED: coax'
  
END 
