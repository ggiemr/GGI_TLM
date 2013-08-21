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
! SUBROUTINE read_output_volume_list
!
! NAME
!     read_output_volume_list
!
! DESCRIPTION
!     read output volume list
!
! Example packet:
!
!Output_volume_list
!2    		! number of output volumes
!1		! OUTPUT NUMBER
!1		! volume number for output
!Ex
!output_time_information
!0.0	! first output time
!1e-6	! last output time
!1e-7	! output time interval
!2		! OUTPUT NUMBER
!3		! volume number for output
!Hz
!output_timestep_information
!0	! first output timestep
!200	! last output timestep
!10	! output timestep interval
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 11/02/2013 CJS
!
!
SUBROUTINE read_output_volume_list

USE TLM_general
USE file_information
USE geometry
USE TLM_output
USE cell_parameters

IMPLICIT NONE

! local variables

integer	:: output_number
integer	:: read_number

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_output_volume_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_output_volumes
  
  CALL write_line_integer('number of output volumes',n_output_volumes,0,output_to_screen_flag)
  
  if ( allocated( output_volumes ) ) GOTO 9000
  
  ALLOCATE ( output_volumes(1:n_output_volumes) )

  do output_number=1,n_output_volumes
  
    CALL write_line_integer('Reading output number',output_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.output_number) goto 9010
    
    read(input_file_unit,*,err=9005)output_volumes(output_number)%volume_number
      	
    CALL read_field_component(input_file_unit,output_volumes(output_number)%field_component)
      	
    CALL read_output_time_information(input_file_unit,	&
                                      output_volumes(output_number)%specified_timestep_information,	&
                                      output_volumes(output_number)%first_timestep,	&
                                      output_volumes(output_number)%last_timestep,	&
                                      output_volumes(output_number)%timestep_interval,	&
                                      output_volumes(output_number)%specified_time_information,	&
                                      output_volumes(output_number)%first_time,	&
                                      output_volumes(output_number)%last_time,	&			      
                                      output_volumes(output_number)%time_interval )
    
  end do ! next output volume

  CALL write_line('FINISHED: read_output_volume_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating output_volumes:',0,.TRUE.)
     CALL write_line('output_volumes already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading output volume list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading output volume list packet',0,.TRUE.)
     CALL write_line('output volumes should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_output_volume_list
