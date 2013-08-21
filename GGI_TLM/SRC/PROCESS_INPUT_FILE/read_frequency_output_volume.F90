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
! SUBROUTINE read_frequency_output_volume
!
! NAME
!     read_frequency_output_volume
!
! DESCRIPTION
!     read frequency output volume packet
!
! Example packet:
!
!Frequency_output_volume_list
!1    		! number of Frequency output volumes
!1		! FREQUENCY OUTPUT VOLUME NUMBER
!1     volume number
!1.3E7 frequency for output
!Ex
!
! COMMENTS
!     
!
! HISTORY
!
!     started 1/3/2013 CJS
!
!
SUBROUTINE read_frequency_output_volume

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

! local variables

  integer	:: i
  integer	:: read_number
  
! START  

  CALL write_line('CALLED: read_frequency_output_volume',0,output_to_screen_flag)
  
  read(input_file_unit,*,err=9000)n_frequency_output_volumes
  
  CALL write_line_integer('number of frequency output volumes',n_frequency_output_volumes,0,output_to_screen_flag)
  
  if ( allocated( frequency_output_volume ) ) GOTO 9010
  
  ALLOCATE( frequency_output_volume(1:n_frequency_output_volumes) )
  
  do i=1,n_frequency_output_volumes
  
    CALL write_line_integer('Reading requency_output_volume number',i,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9000)read_number
    if (read_number.ne.i) GOTO 9020
  
    read(input_file_unit,*,err=9000)frequency_output_volume(i)%volume_number
       
    read(input_file_unit,*,err=9000)frequency_output_volume(i)%frequency
      	
    CALL read_field_component(input_file_unit,frequency_output_volume(i)%field_component)
     
  end do ! next frequency output volume    
  
  CALL write_line('FINISHED: read_frequency_output_volume',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading frequency_output_volume packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9010 CALL write_line('Error allocating frequency_output_volume:',0,.TRUE.)
     CALL write_line('frequency_output_volume already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading frequency_output_volume packet',0,.TRUE.)
     CALL write_line('frequency output volumes should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     STOP
  
END SUBROUTINE read_frequency_output_volume
