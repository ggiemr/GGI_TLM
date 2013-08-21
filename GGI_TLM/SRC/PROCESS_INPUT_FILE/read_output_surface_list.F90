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
! SUBROUTINE read_output_surface_list
!
! NAME
!     read_output_surface_list
!
! DESCRIPTION
!     read output surface list
!
! Example packet:
!
!Output_surface_list
!2    		! number of output surfaces
!1		! OUTPUT NUMBER
!1		! surface number for output
!1		! side of surface for output (+1 or -1)
!Ex
!output_time_information
!0.0	! first output time
!1e-6	! last output time
!1e-7	! output time interval
!2		! OUTPUT NUMBER
!3		! surface number for output
!-1		! side of surface for output (+1 or -1)
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
!     started 13/09/2012 CJS
!
!
SUBROUTINE read_output_surface_list

USE TLM_general
USE file_information
USE geometry
USE TLM_output
USE cell_parameters

IMPLICIT NONE

! local variables

integer	:: output_number
integer	:: read_number
integer	:: side_of_surface_for_output

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_output_surface_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_output_surfaces
  
  CALL write_line_integer('number of output surfaces',n_output_surfaces,0,output_to_screen_flag)
  
  if ( allocated( output_surfaces ) ) GOTO 9000
  
  ALLOCATE ( output_surfaces(1:n_output_surfaces) )

  do output_number=1,n_output_surfaces
  
    CALL write_line_integer('Reading output number',output_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.output_number) goto 9010
    
    read(input_file_unit,*,err=9005)output_surfaces(output_number)%surface_number
    read(input_file_unit,*,err=9005)side_of_surface_for_output
    
    if (side_of_surface_for_output.eq.1) then
      output_surfaces(output_number)%output_on_outward_normal=.TRUE.
    else if (side_of_surface_for_output.eq.-1) then
      output_surfaces(output_number)%output_on_outward_normal=.FALSE.
    else 
      GOTO 9020
    end if
      	
    CALL read_field_component_surface_output(input_file_unit,output_surfaces(output_number)%field_component)
      	
    CALL read_output_time_information(input_file_unit,	&
                                      output_surfaces(output_number)%specified_timestep_information,	&
                                      output_surfaces(output_number)%first_timestep,	&
                                      output_surfaces(output_number)%last_timestep,	&
                                      output_surfaces(output_number)%timestep_interval,	&
                                      output_surfaces(output_number)%specified_time_information,	&
                                      output_surfaces(output_number)%first_time,	&
                                      output_surfaces(output_number)%last_time,	&			      
                                      output_surfaces(output_number)%time_interval )
    
  end do ! next output surface

  CALL write_line('FINISHED: read_output_surface_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating output_surfaces:',0,.TRUE.)
     CALL write_line('output_surfaces already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading output surface list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading output surface list packet',0,.TRUE.)
     CALL write_line('output surfaces should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading output surface list packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_output_surface_list
