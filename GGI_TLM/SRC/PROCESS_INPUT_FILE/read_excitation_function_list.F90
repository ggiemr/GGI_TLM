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
! SUBROUTINE read_excitation_function_list
!
! NAME
!     read_excitation_function
!
! DESCRIPTION
!     read excitation function list
!
! Example packet:
!
!Excitation_function_list
!2       ! number of excitations
!1       ! EXCITATION NUMBER
!impulse
!0.0
!2       ! EXCITATION NUMBER
!gaussian
!1e-9 3e-9
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE read_excitation_function_list

USE TLM_general
USE file_information
USE geometry
USE TLM_excitation

IMPLICIT NONE

! local variables

integer	:: excitation_number
integer	:: read_number
integer	:: n_params
integer	:: i

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_excitation_function_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_excitation_functions
  
  CALL write_line_integer('number of excitation functions',n_excitation_functions,0,output_to_screen_flag)
  
  if ( allocated( excitation_functions ) ) GOTO 9000
  
  allocate ( excitation_functions(1:n_excitation_functions) )

  do excitation_number=1,n_excitation_functions
  
    CALL write_line_integer('Reading excitation number',excitation_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.excitation_number) goto 9010
 
! read excitation type string
    read(input_file_unit,'(A)',end=9010)input_line
   
    CALL write_line( '...Reading excitation_function_type:'//trim(input_line),0,.TRUE. )
    
! convert text to lower case
    CALL convert_to_lower_case(input_line,256)
   
    if (input_line.eq.'impulse') then
    
      n_params=1
      excitation_functions(excitation_number)%type=excitation_function_type_impulse
   
    else if (input_line.eq.'gaussian_step') then
    
      n_params=3
      excitation_functions(excitation_number)%type=excitation_function_type_gaussian_step
   
    else if (input_line.eq.'step') then
    
      n_params=2
      excitation_functions(excitation_number)%type=excitation_function_type_step
   
    else if (input_line.eq.'gaussian') then
    
      n_params=3
      excitation_functions(excitation_number)%type=excitation_function_type_gaussian
   
    else if (input_line.eq.'sinusoid') then
    
      n_params=3
      excitation_functions(excitation_number)%type=excitation_function_type_sinusoid
   
    else if (input_line.eq.'gaussian_sinusoid') then
    
      n_params=5
      excitation_functions(excitation_number)%type=excitation_function_type_gaussian_sinusoid
   
    else if (input_line.eq.'gaussian_step_sinusoid') then
    
      n_params=5
      excitation_functions(excitation_number)%type=excitation_function_type_gaussian_step_sinusoid

    else if (input_line.eq.'double_exponential') then
    
      n_params=3
      excitation_functions(excitation_number)%type=excitation_function_type_double_exponential

    else
    
      goto 9030

    end if
    
! read parameters
    excitation_functions(excitation_number)%n_parameters=n_params
    if (n_params.gt.0) then
      read(input_file_unit,*,end=9020)(excitation_functions(excitation_number)%parameters(i),i=1,n_params)
    end if

  end do ! next excitation function

  CALL write_line('FINISHED: read_excitation_function_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating excitation_functions:',0,.TRUE.)
     CALL write_line('excitation_functions already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading excitation function list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Excitation functions should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9020 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Error reading parameters',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9030 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Unrecognised excitation function:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_excitation_function_list
