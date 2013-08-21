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
! SUBROUTINE read_excitation_surface_list
!
! NAME
!     read_excitation_surface
!
! DESCRIPTION
!     read excitation surface list
!
! Example packet:
!
!Excitation_surface_list
!1	       ! number of output surfaces
!1	       ! EXCITATION surface NUMBER
!1	       ! excitation function number
!2	       ! excitation surface
!-1	       ! side of surface for excitation
!Ex
!soft

! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE read_excitation_surface_list

USE TLM_general
USE file_information
USE geometry
USE TLM_excitation

IMPLICIT NONE

! local variables

integer	:: excitation_number
integer	:: read_number
integer	:: side_of_surface_for_excitation

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_excitation_surface_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_excitation_surfaces
  
  CALL write_line_integer('number of excitation surfaces',n_excitation_surfaces,0,output_to_screen_flag)
  
  if ( allocated( excitation_surfaces ) ) GOTO 9000
  
  allocate ( excitation_surfaces(1:n_excitation_surfaces) )

  do excitation_number=1,n_excitation_surfaces
  
    CALL write_line_integer('Reading excitation number',excitation_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.excitation_number) goto 9010

    read(input_file_unit,*,err=9005)excitation_surfaces(excitation_number)%excitation_function_number

    read(input_file_unit,*,err=9005)excitation_surfaces(excitation_number)%surface_number 
       
    read(input_file_unit,*,err=9005)side_of_surface_for_excitation
    
    if (side_of_surface_for_excitation.eq.1) then
      excitation_surfaces(excitation_number)%excitation_on_outward_normal=.TRUE.
    else if (side_of_surface_for_excitation.eq.-1) then
      excitation_surfaces(excitation_number)%excitation_on_outward_normal=.FALSE.
    else 
      GOTO 9020
    end if
    
    CALL read_field_component(input_file_unit,excitation_surfaces(excitation_number)%field_component)

    read(input_file_unit,'(A)',err=9005)input_line
    
    if (input_line.eq.'soft') then
    
      excitation_surfaces(excitation_number)%source_type=source_type_soft   
    
    else if (input_line.eq.'hard') then
    
      excitation_surfaces(excitation_number)%source_type=source_type_hard
    
    else
    
      GOTO 9030
      
    end if
    
  end do ! next excitation surface

  CALL write_line('FINISHED: read_excitation_surface_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating excitation_surfaces:',0,.TRUE.)
     CALL write_line('excitation_surfaces already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading excitation surface list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading excitation surface list packet',0,.TRUE.)
     CALL write_line('Excitation surfaces should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading excitation surface list packet',0,.TRUE.)
     CALL write_line("Side of surface for excitation should be +1 or -1",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9030 CALL write_line('Error reading excitation surface list packet',0,.TRUE.)
     CALL write_line("Excitation type should be 'hard' or 'soft'",0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_excitation_surface_list
