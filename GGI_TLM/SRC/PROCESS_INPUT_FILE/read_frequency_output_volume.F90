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
! Alternative for simple xyz output
!
!Frequency_output_volume_list_xyz
!1    		! number of Frequency output volumes
!output_every 5
!1		! FREQUENCY OUTPUT VOLUME NUMBER
!1     volume number
!1.3E7 frequency for output
!Ex
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 1/3/2013 CJS
!
!
SUBROUTINE read_frequency_output_volume(format_type)

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

integer :: format_type

! local variables

  integer	:: i
  integer	:: read_number
  
  character(LEN=80) :: ipline_local
  character(LEN=80) :: ipline_local2
  
! START  

  CALL write_line('CALLED: read_frequency_output_volume',0,output_to_screen_flag)

  read(input_file_unit,*,err=9000,end=9000)n_frequency_output_volumes
  
  CALL write_line_integer('number of frequency output volumes',n_frequency_output_volumes,0,output_to_screen_flag)
    
  if (format_type.EQ.frequency_output_volume_format_normal) then
  
    frequency_output_volume_format=frequency_output_volume_format_normal
    frequency_output_volume_output_every=1
  	  
  else if (format_type.EQ.frequency_output_volume_format_xyz_field) then
  
    frequency_output_volume_format=frequency_output_volume_format_xyz_field

! read optional 'output_every' line    

    read(input_file_unit,'(A80)',err=9000,end=9000)ipline_local
    
    write(*,*)'Read line:',ipline_local
    
    CALL convert_to_lower_case(ipline_local,80)
    
    if (ipline_local(1:12).EQ.'output_every') then
      ipline_local2=ipline_local(13:LEN(ipline_local))
      read(ipline_local2,*,err=9030,end=9030)frequency_output_volume_output_every
      write(*,*)'Output every=',frequency_output_volume_output_every
    else
      frequency_output_volume_output_every=1
      backspace(input_file_unit)
    end if
    
  else
    
    write(*,*)'Unknown frequency_output_volume_format:',format_type  
    STOP 1
    
  end if
  
  if ( allocated( frequency_output_volume ) ) GOTO 9010
  
  ALLOCATE( frequency_output_volume(1:n_frequency_output_volumes) )
  
  do i=1,n_frequency_output_volumes
  
    CALL write_line_integer('Reading frequency_output_volume number',i,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9000,end=9000)read_number
    if (read_number.ne.i) GOTO 9020
  
    read(input_file_unit,*,err=9000,end=9000)frequency_output_volume(i)%volume_number
       
    read(input_file_unit,*,err=9000,end=9000)frequency_output_volume(i)%frequency
      	
    CALL read_field_component(input_file_unit,frequency_output_volume(i)%field_component)
    
    CALL read_output_time_information(input_file_unit,	&
                                      frequency_output_volume(i)%specified_timestep_information,	&
                                      frequency_output_volume(i)%first_timestep,	&
                                      frequency_output_volume(i)%last_timestep,	&
                                      frequency_output_volume(i)%timestep_interval,	&
                                      frequency_output_volume(i)%specified_time_information,	&
                                      frequency_output_volume(i)%first_time,	&				      
                                      frequency_output_volume(i)%last_time,	&			      
                                      frequency_output_volume(i)%time_interval )
     
  end do ! next frequency output volume    
  
  CALL write_line('FINISHED: read_frequency_output_volume',0,output_to_screen_flag)
    
  RETURN
    
9000 CALL write_line('Error reading frequency_output_volume packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
  
9010 CALL write_line('Error allocating frequency_output_volume:',0,.TRUE.)
     CALL write_line('frequency_output_volume already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9020 CALL write_line('Error reading frequency_output_volume packet',0,.TRUE.)
     CALL write_line('frequency output volumes should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9030 CALL write_line('Error reading output_every line',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
  
END SUBROUTINE read_frequency_output_volume
