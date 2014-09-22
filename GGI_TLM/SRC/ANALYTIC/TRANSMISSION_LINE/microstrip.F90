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
! Name microstrip
!     
!
! Description
!     calculate the impedance and effective permittivity of a mictrostrp
!     given the geometric and dielectric parameters
!     Based on formaule from IPC-2141A (2004) "Design Guide for High-Speed Controlled Impedance Circuit Boards"
!
! Comments:
!      checked against the web based microstrip calculator at http://www.eeweb.com/toolbox/microstrip-impedance
!
! History
!
!     started 19/9/2014 CJS
!

PROGRAM microstrip

USE constants

IMPLICIT NONE

! local_variables

  real w,t,h,er
  real er_eff
  real w_eff
  real e
  real Z_microstrip
  real C_microstrip
  real L_microstrip
  
  real t1,t2
 
! function_types

! START

!  write(*,*)'CALLED: microstrip'

  write(*,*)'Enter microstrip width (mm)'
  read(*,*)w
  write(*,*)'Enter microstrip thickness (mm)'
  read(*,*)t
  write(*,*)'Enter dielectric thickness (mm)'
  read(*,*)h
  write(*,*)'Enter dielectric relative permittivity '
  read(*,*)er
  
  if (w/h.lt.1.0) then
  
    er_eff=(er+1.0)/2.0+(er-1.0)/2.0*	&
           ( sqrt(w/(w+12*h))+0.04*((1.0-w/h)**2) )
	   
  else
  
    er_eff=(er+1.0)/2.0+(er-1.0)/2.0*	&
           ( sqrt(w/(w+12*h)) )
  
  end if
  
  write(*,*)
  write(*,*)'Effective permittivity=',er_eff
  
  e=exp(1.0)
  
  w_eff=w+(t/pi)*log( 4.0*e/( (t/h)**2+(t/(w*pi+1.1*t*pi))**2) )*	&
                (er_eff+1.0)/(2.0*er_eff)

  write(*,*)'Effective width=',w_eff,'mm'
  
  t1=      4.0*( (14.0*er_eff+8.0)/(11.0*er_eff) )*(h/w_eff)
  t2=sqrt(16.0*( (14.0*er_eff+8.0)/(11.0*er_eff) )**2*(h/w_eff)**2+( (er_eff+1.0)/(2.0*er_eff) )*pi**2 )
  
  Z_microstrip=( z0/(2.0*pi*sqrt(2.0*(er_eff+1.0))) )*	&
        log( 1.0+4.0*(h/w_eff)*( t1 + t2 ) )
  
  write(*,*)'Impedance=',Z_microstrip,'ohms'
  
  C_microstrip=sqrt(er_eff)/(c0*Z_microstrip)
  
  write(*,*)'Capacitance per unit length=',C_microstrip,'F/m'
  
  L_microstrip=Z_microstrip**2*C_microstrip
  
  write(*,*)'Inductance per unit length =',L_microstrip,'H/m'
  
  write(*,*)'Wave velocity =',1.0/sqrt(L_microstrip*C_microstrip),'m/s'
  

!  write(*,*)'FINISHED: microstrip'
  
END 

  
  
  
  
