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
! SUBROUTINE read_excitation_mode_list
!
! NAME
!     read_excitation_mode_list
!
! DESCRIPTION
!     read surface list packet
!
! Example packet:
!
!
!Excitation_mode_list
!1    		! number of excitation modes
!1    		! EXCITATION MODE NUMBER
!1		! excitation function number
!1		! surface number
!1		! side of surface for excitation
!Ex
!soft
!waveguide_original_TLM.frequency_surface_field.fout
!1                ! x column 
!2                ! y column 
!3                ! z column 
!6                ! data column to use

! COMMENTS
!     
!
! HISTORY
!
!     started 10/01/2013 CJS
!
!
SUBROUTINE read_excitation_mode_list

USE TLM_general
USE file_information
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

integer	:: mode_number
integer :: read_number
integer :: side_of_surface_for_excitation

logical	:: file_exists

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_excitation_mode_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_excitation_modes
  
  CALL write_line_integer('number of excitation modes',n_excitation_modes,0,output_to_screen_flag)
  
  if ( allocated( excitation_mode_list ) ) GOTO 9000
  
  allocate ( excitation_mode_list(1:n_excitation_modes) )

  do mode_number=1,n_excitation_modes
  
    CALL write_line_integer('Reading mode number',mode_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.mode_number) goto 9010
 
    read(input_file_unit,*,err=9005)excitation_mode_list(mode_number)%excitation_function_number
    read(input_file_unit,*,err=9005)excitation_mode_list(mode_number)%surface_number
       
    read(input_file_unit,*,err=9005)side_of_surface_for_excitation
    
    if (side_of_surface_for_excitation.eq.1) then
      excitation_mode_list(mode_number)%excitation_on_outward_normal=.TRUE.
    else if (side_of_surface_for_excitation.eq.-1) then
      excitation_mode_list(mode_number)%excitation_on_outward_normal=.FALSE.
    else 
      GOTO 9020
    end if
    
    CALL read_field_component(input_file_unit,excitation_mode_list(mode_number)%field_component)

    read(input_file_unit,'(A)',err=9005)input_line
    
    if (input_line.eq.'soft') then
    
      excitation_mode_list(mode_number)%source_type=source_type_soft   
    
    else if (input_line.eq.'hard') then
    
      excitation_mode_list(mode_number)%source_type=source_type_hard
    
    else
    
      GOTO 9030
      
    end if
 
    read(input_file_unit,'(A)',err=9005)excitation_mode_list(mode_number)%mode_file_name 

    CALL write_line('Checking the existance of file:',0,.TRUE.)
    CALL write_line(trim(excitation_mode_list(mode_number)%mode_file_name),0,.TRUE.)
      
    inquire(file=excitation_mode_list(mode_number)%mode_file_name,EXIST=file_exists)
      
    if (.NOT.file_exists) then
! error - no mode field file exists
      goto 9040
    end if
    
    read(input_file_unit,*,err=9005)excitation_mode_list(mode_number)%xcol
    read(input_file_unit,*,err=9005)excitation_mode_list(mode_number)%ycol
    read(input_file_unit,*,err=9005)excitation_mode_list(mode_number)%zcol
    read(input_file_unit,*,err=9005)excitation_mode_list(mode_number)%mode_col
  
  end do ! next mode in Excitation_mode_list

  CALL write_line('FINISHED: read_excitation_mode_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating Excitation_mode_list:',0,.TRUE.)
     CALL write_line('Excitation_mode_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading excitation_mode_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9010 CALL write_line('Error reading excitation_mode_list packet from input file:',0,.TRUE.)
     CALL write_line('Excitation modes should be numbered in order at the moment...',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading excitation_mode_list packet',0,.TRUE.)
     CALL write_line("Side of surface for excitation should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9030 CALL write_line('Error reading excitation_mode_list packet',0,.TRUE.)
     CALL write_line("Excitation type should be 'hard' or 'soft'",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9040 CALL write_line('Error reading excitation_mode_list packet',0,.TRUE.)
     CALL write_line('Mode field file does not exist. Filename:',0,.TRUE.)
     CALL write_line(trim(excitation_mode_list(mode_number)%mode_file_name),0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_excitation_mode_list
