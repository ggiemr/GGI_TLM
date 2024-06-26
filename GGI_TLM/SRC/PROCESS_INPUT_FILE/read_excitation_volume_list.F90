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
! SUBROUTINE read_excitation_volume_list
!
! NAME
!     read_excitation_volume
!
! DESCRIPTION
!     read excitation volume list
!
! Example packet:
!
!Excitation_volume_list
!1	       ! number of output volumes
!1	       ! EXCITATION VOLUME NUMBER
!1	       ! excitation function number
!2	       ! excitation volume
!Ex
!soft

! COMMENTS
!     
!
! HISTORY
!
!     started 10/05/2017 CJS based on read_excitation_surface_list.F90
!
!
SUBROUTINE read_excitation_volume_list

USE TLM_general
USE file_information
USE geometry
USE TLM_excitation

IMPLICIT NONE

! local variables

integer	:: excitation_number
integer	:: read_number
integer	:: side_of_volume_for_excitation

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_excitation_volume_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005,end=9005)n_excitation_volumes
  
  CALL write_line_integer('number of excitation volumes',n_excitation_volumes,0,output_to_screen_flag)
  
  if ( allocated( excitation_volumes ) ) GOTO 9000
  
  allocate ( excitation_volumes(1:n_excitation_volumes) )

  do excitation_number=1,n_excitation_volumes
  
    CALL write_line_integer('Reading excitation number',excitation_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005,end=9005)read_number
    if (read_number.ne.excitation_number) goto 9010

    read(input_file_unit,*,err=9005,end=9005)excitation_volumes(excitation_number)%excitation_function_number

    read(input_file_unit,*,err=9005,end=9005)excitation_volumes(excitation_number)%volume_number 
           
    CALL read_field_component(input_file_unit,excitation_volumes(excitation_number)%field_component)

    read(input_file_unit,'(A)',err=9005,end=9005)input_line
    
    if (input_line.eq.'soft') then
    
      excitation_volumes(excitation_number)%source_type=source_type_soft   
    
    else if (input_line.eq.'hard') then
    
      excitation_volumes(excitation_number)%source_type=source_type_hard
    
    else
    
      GOTO 9030
      
    end if
    
  end do ! next excitation volume

  CALL write_line('FINISHED: read_excitation_volume_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating excitation_volumes:',0,.TRUE.)
     CALL write_line('excitation_volumes already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
  
9005 CALL write_line('Error reading excitation volume list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9010 CALL write_line('Error reading excitation volume list packet',0,.TRUE.)
     CALL write_line('Excitation volumes should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9020 CALL write_line('Error reading excitation volume list packet',0,.TRUE.)
     CALL write_line("Side of volume for excitation should be 0, +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9030 CALL write_line('Error reading excitation volume list packet',0,.TRUE.)
     CALL write_line("Excitation type should be 'hard' or 'soft'",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
  
END SUBROUTINE read_excitation_volume_list
