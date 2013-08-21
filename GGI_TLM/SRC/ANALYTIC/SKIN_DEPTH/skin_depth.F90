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
  program skin_depth
  
! calculate the skin depth in conductors
    
USE constants   

! local variables

  real*8 	:: sigma,mur,f
  
  real*8	:: delta
      
! START

  write(*,*)'Enter conductivity (S/m)'
  read(*,*)sigma
  write(*,*)'Enter relative permeability'
  read(*,*)mur
  write(*,*)'Enter frequency'
  read(*,*)f
  
  delta=sqrt(2d0/(2d0*pi*f*mu0*mur*sigma))
  
  write(*,*)'Skin depth=',delta,' m'
   
 
  end
