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
! SUBROUTINE read_frequency_domain_power_surface
!
! NAME
!     read_frequency_domain_power_surface
!
! DESCRIPTION
!     read frequency domain power surface list packet
!
! Example packet:
!
!Frequency_domain_power_surface_list 	
!1               ! number of frequency domain power surfaces
!1		 ! FREQUENCY DOMAIN POWER SURFACE NUMBER 
!1		 ! geometric surface number
!1		 ! side of surface for output
!100e6 25e9 100  ! fmin fmax n_frequencies for power output
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/2/2013 CJS
!
!
SUBROUTINE read_frequency_domain_power_surface

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

! local variables

  integer	:: i
  integer	:: read_number
  integer	:: side_of_surface_for_output
  
! START  

  CALL write_line('CALLED: read_frequency_domain_power_surface',0,output_to_screen_flag)
  
  read(input_file_unit,*,err=9000)n_frequency_domain_power_surfaces
  
  CALL write_line_integer('number of frequency output surfaces',n_frequency_domain_power_surfaces,0,output_to_screen_flag)
  
  if ( allocated( frequency_domain_power_surface ) ) GOTO 9010
  
  ALLOCATE( frequency_domain_power_surface(1:n_frequency_domain_power_surfaces) )
  
  do i=1,n_frequency_domain_power_surfaces
  
    CALL write_line_integer('Reading frequency_domain_surface number',i,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9000)read_number
    if (read_number.ne.i) GOTO 9020
  
    read(input_file_unit,*,err=9000)frequency_domain_power_surface(i)%surface_number
    
    read(input_file_unit,*,err=9000)side_of_surface_for_output
    
    if (side_of_surface_for_output.eq.1) then
      frequency_domain_power_surface(i)%output_on_outward_normal=.TRUE.
    else if (side_of_surface_for_output.eq.-1) then
      frequency_domain_power_surface(i)%output_on_outward_normal=.FALSE.
    else 
      GOTO 9030
    end if
    
    read(input_file_unit,*,err=9000)	frequency_domain_power_surface(i)%fmin,	&
    					frequency_domain_power_surface(i)%fmax,	&
    					frequency_domain_power_surface(i)%n_frequencies
					
    if (i.ne.1) then
      if (frequency_domain_power_surface(i)%n_frequencies.NE.frequency_domain_power_surface(1)%n_frequencies) then
        GOTO 9040
      end if
    end if 
    
    frequency_domain_power_surface(i)%fstep=							&
       ( (frequency_domain_power_surface(i)%fmax-frequency_domain_power_surface(i)%fmin)	&
                         /dble(frequency_domain_power_surface(i)%n_frequencies-1) )
     
  end do ! next frequency output surface    
  
  CALL write_line('FINISHED: read_frequency_domain_power_surface',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading frequency_domain_power_surface packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9010 CALL write_line('Error allocating frequency_domain_power_surface:',0,.TRUE.)
     CALL write_line('frequency_domain_power_surface already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading frequency_domain_power_surface packet',0,.TRUE.)
     CALL write_line('frequency output surfaces should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9030 CALL write_line('Error reading frequency_domain_power_surface packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9040 CALL write_line('Error reading frequency_domain_power_surface packet',0,.TRUE.)
     CALL write_line("All surfaces should have the same number of frequencies at the moment...",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_frequency_domain_power_surface
