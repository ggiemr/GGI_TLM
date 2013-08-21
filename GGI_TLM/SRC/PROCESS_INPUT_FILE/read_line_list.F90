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
! SUBROUTINE read_line_list
!
! NAME
!     read_line_list
!
! DESCRIPTION
!     read line list packet
!
! Example packet:
!
!line_list
!2   Number of lines (integer)
!1       line_number (integer)
!straight_line
!2.0      line parameters (n*real)
!0.0 0.0 45.0
!-1.5 -0.5 0.0
!2       line_number (integer)
!straight_line2
!0.0 0.0 0.0     1.0 2.0 3.0      line parameters (n*real)
!0.0 0.0 0.0
!0.0 0.0 0.0
!
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
SUBROUTINE read_line_list

USE TLM_general
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: line_number
integer :: read_number

integer	:: n_params
integer :: i

character*256	:: input_line
character*256	:: tri_line_filename

logical	:: file_exists

! START  

  CALL write_line('CALLED: Read_line_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_lines
  
  CALL write_line_integer('number of lines',n_lines,0,output_to_screen_flag)
  
  if ( allocated( problem_lines ) ) GOTO 9000
  
  allocate ( problem_lines(1:n_lines) )

  do line_number=1,n_lines
  
    CALL write_line_integer('Reading line number',line_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.line_number) goto 9010
    
    problem_lines(line_number)%line_number=read_number
 
! read line type string
    read(input_file_unit,'(A)',end=9010)input_line
   
    CALL write_line( '...STARTED reading line_type:'//trim(input_line),0,.TRUE. )
    
! convert text to lower case
    CALL convert_to_lower_case(input_line,256)
   
    if (input_line.eq.'straight_line') then
      n_params=1
      problem_lines(line_number)%line_type=line_type_straight_line

    else if (input_line.eq.'straight_line2') then
      n_params=6
      problem_lines(line_number)%line_type=line_type_straight_line2

    else if (input_line.eq.'arc') then
      n_params=3
      problem_lines(line_number)%line_type=line_type_arc
	      
    else
    
      goto 9030

    end if
    
! read parameters
    problem_lines(line_number)%n_line_parameters=n_params
    if (n_params.gt.0) then
      read(input_file_unit,*,end=9020)    &
      (problem_lines(line_number)%line_parameters(i),i=1,n_params)
    end if
        
! read transformation
    problem_lines(line_number)%trans%trans_type='euler'
    problem_lines(line_number)%trans%trans_number=1
    read(input_file_unit,*,end=9040)    &
    (problem_lines(line_number)%trans%parameters(i),i=1,6)
    
! convert euler angles to radians
    problem_lines(line_number)%trans%parameters(1:3)=pi/180d0	&
       *problem_lines(line_number)%trans%parameters(1:3)
   
    CALL write_line( '...FINISHED reading line_type:'//trim(input_line),0,.TRUE. )
  
  end do ! next line in line_list


  CALL write_line('FINISHED: Read_line_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating problem_lines:',0,.TRUE.)
     CALL write_line('problem_lines already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading line_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading line_list_packet_data',0,.TRUE.)
     CALL write_line('lines should be numbered in order at the moment...',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading line_list_packet_data',0,.TRUE.)
     CALL write_line_integer('Number of paramters expected=',n_params,0,.TRUE.)
     STOP
  
9030 CALL write_line('Error reading line_list_packet_data',0,.TRUE.)
     CALL write_line('Unknown line type:'//trim(input_line),0,.TRUE.)
     STOP
  
9040 CALL write_line('Error reading line_list_packet_data',0,.TRUE.)
     CALL write_line('Error reading transformation data',0,.TRUE.)
     STOP
  
9050 CALL write_line('Error in read_line_list_packet_data',0,.TRUE.)
     CALL write_line('Triangulated line file not found',0,.TRUE.)
     CALL write_line(trim(problem_lines(line_number)%filename),0,.TRUE.)
     STOP
  
END SUBROUTINE read_line_list
