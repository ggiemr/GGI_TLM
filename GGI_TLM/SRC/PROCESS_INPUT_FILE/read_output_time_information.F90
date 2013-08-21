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
! SUBROUTINE read_output_time_information
!
! NAME
!     read_output_time_information
!
! DESCRIPTION
!     
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 13/09/2012 CJS
!
!
SUBROUTINE read_output_time_information(file_unit,	&
                                        specified_timestep_information,	   &
                                        first_timestep,	   &
                                        last_timestep,	   &
                                        timestep_interval,   &
                                        specified_time_information,  &
                                        first_time,  &					 
                                        last_time,   &				 
                                        time_interval )

IMPLICIT NONE

  integer		:: file_unit

  logical		:: specified_timestep_information
  integer		:: first_timestep
  integer		:: last_timestep
  integer		:: timestep_interval
  
  logical		:: specified_time_information
  real*8		:: first_time
  real*8		:: last_time
  real*8		:: time_interval

! local variables

character*256	:: input_line

! START  
    read(file_unit,'(A)',err=9000)input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,256)
    
    if      (input_line.eq.'output_timestep_information') then   
    
      specified_timestep_information=.TRUE.
      specified_time_information=.FALSE.
      read(file_unit,*)first_timestep
      read(file_unit,*)last_timestep
      read(file_unit,*)timestep_interval
      
    else if (input_line.eq.'output_time_information') then   
    
      specified_time_information=.TRUE.
      specified_timestep_information=.FALSE.
      read(file_unit,*)first_time
      read(file_unit,*)last_time
      read(file_unit,*)time_interval

    else

! set the default output time information and return
      specified_timestep_information=.TRUE.
      specified_time_information=.FALSE.
      first_timestep=0
      last_timestep=1000000000
      timestep_interval=1
      backspace(file_unit)
      
!      GOTO 9010
      
    end if

    
  RETURN
  
9000 CALL write_line('Error reading output time information from file:',0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
     
9010 CALL write_line('Error reading output time information',0,.TRUE.)
     CALL write_line("Expecting either 'output_timestep_information' or 'output_time_information'",0,.TRUE.)
     CALL write_error_line(file_unit)
     STOP
  
END SUBROUTINE read_output_time_information
