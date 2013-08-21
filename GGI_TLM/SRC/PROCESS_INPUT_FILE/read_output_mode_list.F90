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
! SUBROUTINE read_output_mode_list
!
! NAME
!     read_output_mode_list
!
! DESCRIPTION
!     read output mode list packet
!
! Example packet:
!
!Output_mode_list
!1    		! number of output modes
!1    		! OUTPUT MODE NUMBER
!1		! surface number
!1		! side of surface for output
!Ex
!waveguide_original_TLM.frequency_surface_field.fout
!1                ! x column 
!2                ! y column 
!3                ! z column 
!6                ! data column to use
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/01/2013 CJS
!
!
SUBROUTINE read_output_mode_list

USE TLM_general
USE file_information
USE TLM_output
USE constants

IMPLICIT NONE

! local variables

integer	:: mode_number
integer :: read_number
integer :: side_of_surface_for_output

logical	:: file_exists

! START  

  CALL write_line('CALLED: read_output_mode_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_output_modes
  
  CALL write_line_integer('number of output modes',n_output_modes,0,output_to_screen_flag)
  
  if ( allocated( output_mode_list ) ) GOTO 9000
  
  allocate ( output_mode_list(1:n_output_modes) )

  do mode_number=1,n_output_modes
  
    CALL write_line_integer('Reading mode number',mode_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.mode_number) goto 9010
 
    read(input_file_unit,*,err=9005)output_mode_list(mode_number)%surface_number
       
    read(input_file_unit,*,err=9005)side_of_surface_for_output
    
    if (side_of_surface_for_output.eq.1) then
      output_mode_list(mode_number)%output_on_outward_normal=.TRUE.
    else if (side_of_surface_for_output.eq.-1) then
      output_mode_list(mode_number)%output_on_outward_normal=.FALSE.
    else 
      GOTO 9020
    end if
     
    CALL read_field_component(input_file_unit,output_mode_list(mode_number)%field_component)

    read(input_file_unit,'(A)',err=9005)output_mode_list(mode_number)%mode_file_name 

    CALL write_line('Checking the existance of file:',0,.TRUE.)
    CALL write_line(trim(output_mode_list(mode_number)%mode_file_name),0,.TRUE.)
      
    inquire(file=output_mode_list(mode_number)%mode_file_name,EXIST=file_exists)
      
    if (.NOT.file_exists) then
! error - no mode field file exists
      goto 9030
    end if
    
    read(input_file_unit,*,err=9005)output_mode_list(mode_number)%xcol
    read(input_file_unit,*,err=9005)output_mode_list(mode_number)%ycol
    read(input_file_unit,*,err=9005)output_mode_list(mode_number)%zcol
    read(input_file_unit,*,err=9005)output_mode_list(mode_number)%mode_col
    
  end do ! next mode in output_mode_list

  CALL write_line('FINISHED: read_output_mode_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating output_mode_list:',0,.TRUE.)
     CALL write_line('output_mode_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading output_mode_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9010 CALL write_line('Error reading output_mode_list packet from input file:',0,.TRUE.)
     CALL write_line('output modes should be numbered in order at the moment...',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading output_mode_list packet',0,.TRUE.)
     CALL write_line("Side of surface for output should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9030 CALL write_line('Error reading output_mode_list packet',0,.TRUE.)
     CALL write_line('Mode field file does not exist. Filename:',0,.TRUE.)
     CALL write_line(trim(output_mode_list(mode_number)%mode_file_name),0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  

END SUBROUTINE read_output_mode_list
