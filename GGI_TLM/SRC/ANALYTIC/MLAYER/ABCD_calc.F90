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
!  SUBROUTINE ABCD_calc_material
!  SUBROUTINE ABCD_calc_thin_layer

! NAME
!     SUBROUTINE  ABCD_calc_material
!
! DESCRIPTION
! 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
  SUBROUTINE ABCD_calc_material(layer)

USE constants
USE filter_types
USE filter_operators
USE filter_functions

USE Mlayer_file_module
USE Mlayer_module

! variables passed to subroutine

  integer layer

! local variables

  real*8 d
  complex*16 epsr,mur
  complex*16 beta,Z
  
  complex*16 j_beta_d,p_j_beta_d,m_j_beta_d
  real*8     re,im

! START

  d=layer_thickness(layer)

  epsr=evaluate_Sfilter_frequency_response(material_list(layer)%eps_S,f)  &
       -j*material_list(layer)%sigma_e/(w*eps0) 
       
  mur=evaluate_Sfilter_frequency_response(material_list(layer)%mu_S,f)  &
      -j*material_list(layer)%sigma_m/(w*eps0) 
           
  beta=w*sqrt(mu0*mur*eps0*epsr)/cos(angle_rad)
  
  if (aimag(beta).gt.0d0) beta=beta*(-1d0,0d0)
  
  if (polarisation.eq.TE) then
    Z=sqrt(mu0*mur/(eps0*epsr))*cos(angle_rad)
  else
    Z=sqrt(mu0*mur/(eps0*epsr))/cos(angle_rad)
  end if
  
  p_j_beta_d= j*beta*d
  m_j_beta_d=-j*beta*d
  
  re=dble(p_j_beta_d)
  im=imag(p_j_beta_d)
  if (re.gt.20d0) then
    re=20d0
    p_j_beta_d=dcmplx(re,im)
    overflow_error=.TRUE.
  end if
  
  re=dble(m_j_beta_d)
  im=imag(m_j_beta_d)
  if (re.gt.20d0) then
    re=20d0
    m_j_beta_d=dcmplx(re,im)
    overflow_error=.TRUE.
  end if
  
  ABCD(1,1,layer)=0.5d0*  ( exp(p_j_beta_d)+exp(m_j_beta_d) )
  ABCD(1,2,layer)=0.5d0*Z*( exp(p_j_beta_d)-exp(m_j_beta_d) )
  ABCD(2,1,layer)=0.5d0*  ( exp(p_j_beta_d)-exp(m_j_beta_d) )/Z
  ABCD(2,2,layer)=0.5d0*  ( exp(p_j_beta_d)+exp(m_j_beta_d) )
  
         
  return

  END SUBROUTINE ABCD_calc_material
  
  
! NAME
!     SUBROUTINE  ABCD_calc_thin_layer
!
! DESCRIPTION
! 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
  SUBROUTINE ABCD_calc_thin_layer(layer)

USE constants
USE filter_types
USE filter_operators
USE filter_functions

USE Mlayer_file_module
USE Mlayer_module

! variables passed to subroutine

  integer layer

! local variables

  complex*16 Z11,Z12,Z21,Z22

! START

  Z11=evaluate_Sfilter_frequency_response(thin_layer_list(layer)%Z11_S,f)
  Z12=evaluate_Sfilter_frequency_response(thin_layer_list(layer)%Z12_S,f)
  Z21=evaluate_Sfilter_frequency_response(thin_layer_list(layer)%Z21_S,f)
  Z22=evaluate_Sfilter_frequency_response(thin_layer_list(layer)%Z22_S,f)

  Z11=-Z11
  Z21=-Z21
  
  ABCD(1,1,layer)=Z22/Z12
  ABCD(1,2,layer)=Z21-Z11*Z22/Z12
  ABCD(2,1,layer)=1d0/Z12
  ABCD(2,2,layer)=-Z11/Z12

  return

  END SUBROUTINE ABCD_calc_thin_layer
  
  
