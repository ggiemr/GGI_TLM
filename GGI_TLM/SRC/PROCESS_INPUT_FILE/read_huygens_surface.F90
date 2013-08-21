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
! Example packet:
!
!Huygens_surface
!0     surface number
!-1	       ! side of surface for excitation
!1     excitation function number
!0.0 0.0  wave vector Theta and Phi 
!1.0 0.0  Polarisation theta and Phi
!
!
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
    
    read(input_file_unit,*,err=9000)huygens_surface%Ktheta,huygens_surface%Kphi
    
    read(input_file_unit,*,err=9000)huygens_surface%Etheta,huygens_surface%Ephi

! convert angles to radians
    
  huygens_surface%Ktheta=(pi/180d0)*huygens_surface%Ktheta
  huygens_surface%Kphi  =(pi/180d0)*huygens_surface%Kphi
  
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
