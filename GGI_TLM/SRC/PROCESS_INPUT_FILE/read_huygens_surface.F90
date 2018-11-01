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
! SUBROUTINE read_huygens_surface
!
! NAME
!     read_huygens_surface
!
! DESCRIPTION
!     read huygens surface data
!
! Example packet1:
!
!Huygens_surface
!0       !surface number
!-1	 ! side of surface for excitation
!1       !  excitation function number
!0.0 0.0 !  wave vector Theta and Phi 
!1.0 0.0 !  Polarisation theta and Phi
!
!
! Example packet2:
!
!Huygens_surface
!0       !surface number
!-1	 ! side of surface for excitation
!1       !  excitation function number
!Random  !  wave vector Theta and Phi 
!Random  !  Polarisation theta and Phi
!!
! COMMENTS
!     
!
! HISTORY
!
!     started 29/11/2012 CJS
!
!
SUBROUTINE read_huygens_surface

USE TLM_general
USE file_information
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

  integer side_of_surface_for_excitation
  
  real*8 a,b
  real*8 r,theta,phi,u
  real*8 x,y,z
  character ch

! START  

  CALL write_line('CALLED: read_huygens_surface',0,output_to_screen_flag)
  
    n_huygens_surfaces=1
    
    read(input_file_unit,*,err=9000)huygens_surface%surface_number
    
    if (huygens_surface%surface_number.EQ.0) then
      huygens_surface%outer_surface_flag=.TRUE.
    else
      huygens_surface%outer_surface_flag=.FALSE.
    end if
       
    read(input_file_unit,*,err=9000)side_of_surface_for_excitation
    
    if (side_of_surface_for_excitation.eq.1) then
      huygens_surface%excitation_on_outward_normal=.TRUE.
    else if (side_of_surface_for_excitation.eq.-1) then
      huygens_surface%excitation_on_outward_normal=.FALSE.
    else 
      GOTO 9020
    end if
    
    read(input_file_unit,*,err=9000)huygens_surface%excitation_function_number

! check for random excitation    
    
! Wave direction 
    read(input_file_unit,'(A1)')ch
    
    if ( (ch.eq.'r').OR.(ch.eq.'R') ) then

! set random direction
      CALL random_number(a)
      CALL random_number(b)
      theta=a*2d0*pi
      u    =b*2d0-1d0
    
      x=cos(theta)*sqrt(1d0-u*u)
      y=sin(theta)*sqrt(1d0-u*u)
      z=u
    
      r=sqrt(x*x+y*y+z*z)
      if (x.ne.0d0) then
        phi=atan2(y,x)
      else
        if (y.gt.0d0) then
          phi=pi/2d0
        else if (y.lt.0d0) then
          phi=-pi/2d0
        else
          phi=0d0
        end if
      end if
      if(r.ne.0d0) then
        theta=acos(z/r)
      else
        theta=0d0
      end if
      
      huygens_surface%Ktheta=theta
      huygens_surface%Kphi=phi

    else
! this is not a face junction so go back as if the extra line had not been read
  
      backspace(unit=input_file_unit)
      read(input_file_unit,*,err=9000)huygens_surface%Ktheta,huygens_surface%Kphi
  
! convert angles to radians
      huygens_surface%Ktheta=(pi/180d0)*huygens_surface%Ktheta
      huygens_surface%Kphi  =(pi/180d0)*huygens_surface%Kphi
  
    end if
       
! Wave polarisation 
    read(input_file_unit,'(A1)')ch
    
    if ( (ch.eq.'r').OR.(ch.eq.'R') ) then

! set random polarisation
  
      CALL random_number(a)
! NEW: convert a to an angle between 0 and 2pi then set Etheta, Ephi as cos and sin of this angle
      a=a*2d0*pi
      huygens_surface%Etheta=cos(a)
      huygens_surface%Ephi=sin(a)
      
! OLD: Etheta, Ephi in first quadrant only
!      b=sqrt(1d0-a*a)
!      huygens_surface%Etheta=a
!      huygens_surface%Ephi=b
      
    else
! this is not a face junction so go back as if the extra line had not been read
  
      backspace(unit=input_file_unit)
      read(input_file_unit,*,err=9000)huygens_surface%Etheta,huygens_surface%Ephi
  
    end if
  
  CALL write_line('FINISHED: read_huygens_surface',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading huygens surface packet data',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading huygens surface packet data',0,.TRUE.)
     CALL write_line("Side of surface for excitation should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_huygens_surface
