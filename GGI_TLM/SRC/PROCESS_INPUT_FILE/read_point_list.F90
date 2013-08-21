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
! SUBROUTINE read_point_list
!
! NAME
!     read_point_list
!
! DESCRIPTION
!     read point list packet
!
! Example packet:
!
!Point_list
!1   Number of points (integer)
!1       point_number (integer)
!1.0 1.0 1.0    point coordinates (3*real)
!0.0 0.0 0.0
!0.0 0.0 0.0
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE read_point_list

USE TLM_general
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: point_number
integer :: read_number

integer	:: n_params
integer :: i

character*256	:: input_line

logical	:: file_exists

! START  

  CALL write_line('CALLED: Read_point_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_points
  
  CALL write_line_integer('number of points',n_points,0,output_to_screen_flag)
  
  if ( allocated( problem_points ) ) GOTO 9000
  
  allocate ( problem_points(1:n_points) )

  do point_number=1,n_points
  
    CALL write_line_integer('Reading point number',point_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.point_number) goto 9010
    
! read coordinates
    read(input_file_unit,*,end=9020)problem_points(point_number)%point%x,	&
    				    problem_points(point_number)%point%y,	&
    				    problem_points(point_number)%point%z
    
    problem_points(point_number)%trans%trans_type='euler'
    problem_points(point_number)%trans%trans_number=1
    read(input_file_unit,*,end=9030)    &
    (problem_points(point_number)%trans%parameters(i),i=1,6)
    
! convert euler angles to radians
    problem_points(point_number)%trans%parameters(1:3)=pi/180d0	&
       *problem_points(point_number)%trans%parameters(1:3)
  
  end do ! next point in point_list

  CALL write_line('FINISHED: Read_point_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating problem_points:',0,.TRUE.)
     CALL write_line('problem_points already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading point_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading point_list_packet_data',0,.TRUE.)
     CALL write_line('points should be numbered in order at the moment...',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading point_list_packet_data',0,.TRUE.)
     CALL write_line('Expecting coordinate data (3*real)',0,.TRUE.)
     STOP
  
9030 CALL write_line('Error reading point_list_packet_data',0,.TRUE.)
     CALL write_line('Error reading transformation data',0,.TRUE.)
     STOP
  
END SUBROUTINE read_point_list
