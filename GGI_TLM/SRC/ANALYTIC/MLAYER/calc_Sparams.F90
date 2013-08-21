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
!  SUBROUTINE calc_Sparams
!  SUBROUTINE S_to_Z2

! NAME
!     SUBROUTINE  calc_Sparams
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
  SUBROUTINE calc_Sparams()

USE constants
USE filter_types
USE filter_operators
USE filter_functions

USE Mlayer_file_module
USE Mlayer_module

! local variables

  integer layer

  complex*16 ABCD_m(2,2),ABCD_temp(2,2)
  complex*16 A,B,C,D
  
  complex*16 S_11,S_12,S_21,S_22
  complex*16 Z_11,Z_12,Z_21,Z_22
  
! function variables
  real*8 atan2b
  
  real*8 Z0_inc

! S_TART

! first assemble the ABCD matrix of the layers
  
  do layer=1,m
    if (layer_type(layer).eq.material) then
      call  ABCD_calc_material(layer)
    else
      call  ABCD_calc_thin_layer(layer)    
    end if
  end do

! multiply the ABCD matrices together
  ABCD_m(1,1)=(1d0,0d0)
  ABCD_m(1,2)=(0d0,0d0)
  ABCD_m(2,1)=(0d0,0d0)
  ABCD_m(2,2)=(1d0,0d0)
  
  do layer=1,m
    call cmatmul(ABCD(1,1,layer),2,2,ABCD_m(1,1),2,2,ABCD_temp(1,1),2)
    ABCD_m(:,:)=ABCD_temp(:,:)
  end do

  A=ABCD_m(1,1)
  B=ABCD_m(1,2)
  C=ABCD_m(2,1)
  D=ABCD_m(2,2)
  
! solve for the S_ parameters

  if (polarisation.eq.TE) then
    Z0_inc=Z0*cos(angle_rad)
  else if (polarisation.eq.TM) then
    Z0_inc=Z0/cos(angle_rad)
  else
    write(*,*)'Dont know what the polarisation is...'
    write(*,*) polarisation
  end if

  S_11=-(A-B/Z0_inc-D+C*Z0_inc)/ &
       (A+B/Z0_inc+D+C*Z0_inc)
       
  S_21=A-B/Z0_inc+S_11*(A+B/Z0_inc)
       
  S_12=2d0/(A+B/Z0_inc+C*Z0_inc+D)
  
  S_22=1d0-S_12*(C*Z0_inc+D)
    
! solve for the impedance parameters    
    
  call S_to_Z2(S_11,S_12,S_21,S_22,Z_11,Z_12,Z_21,Z_22,Z0_inc)
   
! write results to file   
   
  write(S11_file_unit,8010)f,angle,real(S_11),imag(S_11),  &
                           abs(S_11),atan2b(imag(S_11),real(S_11))*180D0/pi  &
			   ,-20d0*log10(abs(S_11))
  write(S12_file_unit,8010)f,angle,real(S_12),imag(S_12),  &
                           abs(S_12),atan2b(imag(S_12),real(S_12))*180D0/pi  &
			   ,-20d0*log10(abs(S_12))
  write(S21_file_unit,8010)f,angle,real(S_21),imag(S_21),  &
                           abs(S_21),atan2b(imag(S_21),real(S_21))*180D0/pi  &
			   ,-20d0*log10(abs(S_21))
  write(S22_file_unit,8010)f,angle,real(S_22),imag(S_22),  &
                           abs(S_22),atan2b(imag(S_22),real(S_22))*180D0/pi  &
			   ,-20d0*log10(abs(S_22))

  write(Z11_file_unit,8000)f,real(Z_11),imag(Z_11),  &
                           abs(Z_11),atan2b(imag(Z_11),real(Z_11))*180D0/pi
  write(Z12_file_unit,8000)f,real(Z_12),imag(Z_12),  &
                           abs(Z_12),atan2b(imag(Z_12),real(Z_12))*180D0/pi
  write(Z21_file_unit,8000)f,real(Z_21),imag(Z_21),  &
                           abs(Z_21),atan2b(imag(Z_21),real(Z_21))*180D0/pi
  write(Z22_file_unit,8000)f,real(Z_22),imag(Z_22),  &
                           abs(Z_22),atan2b(imag(Z_22),real(Z_22))*180D0/pi

  write(power_file_unit,8000)f,abs(S_11*conjg(S_11))+abs((S_21)*conjg(s_21)), &
                               abs(S_22*conjg(S_22))+abs((S_12)*conjg(s_12))
			       

8000 format(5E16.8)
8010 format(7E12.4)
  return

  END SUBROUTINE calc_Sparams
  
! NAME
!     SUBROUTINE  S_to_Z2
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
!     
!
!
!
  SUBROUTINE S_to_Z2(S11,S12,S21,S22,Z11,Z12,Z21,Z22,Z0_inc)
  
USE constants
             
  complex*16 S11,S12,S21,S22
  complex*16 Z11,Z12,Z21,Z22
  real*8 Z0_inc
  
  complex*16 IPS(2,2),IMS(2,2),IMSinv(2,2),DET

! START

! calulate impedance parameters Z=Z0 [I-S]^-1 [I+S]	  
  IMS(1,1)=(1D0,0D0)-S11
  IMS(1,2)=	    -S12
  IMS(2,1)=	    -S21
  IMS(2,2)=(1D0,0D0)-S22

  IPS(1,1)=(1D0,0D0)+S11
  IPS(1,2)=	    +S12
  IPS(2,1)=	    +S21
  IPS(2,2)=(1D0,0D0)+S22

  DET=IMS(1,1)*IMS(2,2)-IMS(1,2)*IMS(2,1)

  IMSinv(1,1)= IMS(2,2)/DET
  IMSinv(1,2)=-IMS(1,2)/DET
  IMSinv(2,1)=-IMS(2,1)/DET
  IMSinv(2,2)= IMS(1,1)/DET

  Z11=IMSinv(1,1)*IPS(1,1)+IMSinv(1,2)*IPS(2,1)
  Z12=IMSinv(1,1)*IPS(1,2)+IMSinv(1,2)*IPS(2,2)
  Z21=IMSinv(2,1)*IPS(1,1)+IMSinv(2,2)*IPS(2,1)
  Z22=IMSinv(2,1)*IPS(1,2)+IMSinv(2,2)*IPS(2,2)

  Z11=Z11*z0_inc
  Z12=Z12*z0_inc
  Z21=Z21*z0_inc
  Z22=Z22*z0_inc
  	  
  return

  END SUBROUTINE S_to_Z2
!
!
!

  real*8 function atan2b(a,b)
  
USE constants 
 
  real*8 a,b
  
  if (b.ne.0d0) then
    atan2b=atan2(a,b)
  else
    if (a.gt.0d0) then
      atan2b=pi/2d0
    else if (a.lt.0d0) then
      atan2b=-pi/2d0
    else
      atan2b=0d0
    end if
  end if
  
  return
  end
