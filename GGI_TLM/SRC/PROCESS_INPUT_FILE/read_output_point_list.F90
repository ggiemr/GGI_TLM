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
! SUBROUTINE read_output_point_list
!
! NAME
!     read_output_point_list
!
! DESCRIPTION
!     read output point list
!
! Example packet:
!
!Output_point_list
!2    		! number of output points
!1		! OUTPUT NUMBER
!1		! point number for output
!Ex
!centre
!output_timestep_information
!0	! first output timestep
!200	! last output timestep
!10	! output timestep interval
!2		! OUTPUT NUMBER
!2		! point number for output
!Hy
!zmax
!output_time_information
!0.0	! first output time
!1e-6	! last output time
!1e-7	! output time interval
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE read_output_point_list

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

  CALL write_line('CALLED: read_output_point_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_output_points
  
  CALL write_line_integer('number of output points',n_output_points,0,output_to_screen_flag)
  
  if ( allocated( output_points ) ) GOTO 9000
  
  allocate ( output_points(1:n_output_points) )

  do output_number=1,n_output_points
  
    CALL write_line_integer('Reading output number',output_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.output_number) goto 9010
    
    read(input_file_unit,*,err=9005)output_points(output_number)%point_number
    
    CALL get_point(output_points(output_number)%point_number,	&
                   output_points(output_number)%cell_point%cell,			&
		   output_points(output_number)%point  )
    	
    CALL read_field_component(input_file_unit,output_points(output_number)%field_component)

    CALL read_centre_or_face(input_file_unit,output_points(output_number)%cell_point%point)
    
    CALL read_output_time_information(input_file_unit,	&
                                      output_points(output_number)%specified_timestep_information,	&
                                      output_points(output_number)%first_timestep,	&
                                      output_points(output_number)%last_timestep,	&
                                      output_points(output_number)%timestep_interval,	&
                                      output_points(output_number)%specified_time_information,	&
                                      output_points(output_number)%first_time,	&				      
                                      output_points(output_number)%last_time,	&			      
                                      output_points(output_number)%time_interval )
     
  end do ! next output point

  CALL write_line('FINISHED: read_output_point_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating output_points:',0,.TRUE.)
     CALL write_line('output_points already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading output point list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading output point list packet',0,.TRUE.)
     CALL write_line('output points should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading excitation point list packet',0,.TRUE.)
     CALL write_line("Output type should be 'centre','xmin','xmax','ymin','ymax','zmin','zmax'",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_output_point_list
