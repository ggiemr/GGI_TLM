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
! SUBROUTINE read_volume_list
!
! NAME
!     read_volume_list
!
! DESCRIPTION
!     read volume list packet
!
! Example packet:
!
!	volume_list
!	3   Number of volumes (integer)
!	1       volume_number (integer)
!	sphere 
!	1.0       volume parameters (n*real)
!	0.0 0.0 0.0
!	0.0 0.0 0.0
!	2       volume_number (integer)
!	rectangular_block
!	2.5 2.5 2.5      volume parameters (n*real)
!	0.0 0.0 0.0
!	0.0 0.0 0.0
!	3       volume_number (integer)
!	rectangular_block
!	3.0 3.0 3.0     volume parameters (n*real)
!	0.0 0.0 0.0
!	0.0 0.0 0.0
!
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE read_volume_list

USE TLM_general
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number
integer :: read_number

integer	:: n_params
integer :: i

character*256	:: input_line
character*256	:: tet_volume_filename

logical	:: file_exists

! START  

  CALL write_line('CALLED: Read_volume_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_volumes
  
  CALL write_line_integer('number of volumes',n_volumes,0,output_to_screen_flag)
  
  if ( allocated( problem_volumes ) ) GOTO 9000
  
  ALLOCATE ( problem_volumes(1:n_volumes) )

  do volume_number=1,n_volumes
  
    CALL write_line_integer('Reading volume number',volume_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.volume_number) goto 9010
    
    problem_volumes(volume_number)%volume_number=read_number
 
! read volume type string
    read(input_file_unit,'(A)',end=9010)input_line
   
    CALL write_line( '...STARTED reading volume_type:'//trim(input_line),0,.TRUE. )
    
! convert text to lower case
    CALL convert_to_lower_case(input_line,256)
   
    if (input_line.eq.'rectangular_block2') then
      n_params=6
      problem_volumes(volume_number)%volume_type=volume_type_rectangular_block2

    else if (input_line.eq.'rectangular_block') then
      n_params=3
      problem_volumes(volume_number)%volume_type=volume_type_rectangular_block
     
    else if (input_line.eq.'cylinder') then
      n_params=2
      problem_volumes(volume_number)%volume_type=volume_type_cylinder

    else if (input_line.eq.'sphere') then
      n_params=1
      problem_volumes(volume_number)%volume_type=volume_type_sphere

    else if (input_line.eq.'tet') then
      n_params=12
      problem_volumes(volume_number)%volume_type=volume_type_tet

    else if (input_line.eq.'pyramid_ram') then
      n_params=3
      problem_volumes(volume_number)%volume_type=volume_type_pyramid_ram

    else if (input_line.eq.'pyramid') then
      n_params=2
      problem_volumes(volume_number)%volume_type=volume_type_pyramid

    else if (input_line.eq.'tet_volume_mesh') then
      n_params=1
      problem_volumes(volume_number)%volume_type=volume_type_tet_mesh

! for a tet volume mesh, read the mesh filename and check that the file exists
      
      read(input_file_unit,'(A)',end=9005)tet_volume_filename
      problem_volumes(volume_number)%filename=tet_volume_filename
      
      CALL write_line('Checking the existance of file:',0,.TRUE.)
      CALL write_line(trim(problem_volumes(volume_number)%filename),0,.TRUE.)
      
      inquire(file= problem_volumes(volume_number)%filename,EXIST=file_exists)
      
      if (.NOT.file_exists) then
! error - no tet volume mesh file exists
        goto 9050
      end if
      	      
    else
    
      goto 9030

    end if
    
! read parameters
    problem_volumes(volume_number)%n_volume_parameters=n_params
    if (n_params.gt.0) then
      read(input_file_unit,*,end=9020)    &
      (problem_volumes(volume_number)%volume_parameters(i),i=1,n_params)
    end if
        
! read transformation
    problem_volumes(volume_number)%trans%trans_type='euler'
    problem_volumes(volume_number)%trans%trans_number=1
    read(input_file_unit,*,end=9040)    &
    (problem_volumes(volume_number)%trans%parameters(i),i=1,6)
    
! convert euler angles to radians
    problem_volumes(volume_number)%trans%parameters(1:3)=pi/180d0	&
       *problem_volumes(volume_number)%trans%parameters(1:3)
   
    CALL write_line( '...FINISHED reading volume_type:'//trim(input_line),0,.TRUE. )
  
  end do ! next volume in volume_list


  CALL write_line('FINISHED: Read_volume_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating problem_volumes:',0,.TRUE.)
     CALL write_line('problem_volumes already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading volume_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading volume_list_packet_data',0,.TRUE.)
     CALL write_line('volumes should be numbered in order at the moment...',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading volume_list_packet_data',0,.TRUE.)
     CALL write_line_integer('Number of paramters expected=',n_params,0,.TRUE.)
     STOP
  
9030 CALL write_line('Error reading volume_list_packet_data',0,.TRUE.)
     CALL write_line('Unknown volume type:'//trim(input_line),0,.TRUE.)
     STOP
  
9040 CALL write_line('Error reading volume_list_packet_data',0,.TRUE.)
     CALL write_line('Error reading transformation data',0,.TRUE.)
     STOP
  
9050 CALL write_line('Error in read_volume_list_packet_data',0,.TRUE.)
     CALL write_line('Tet volume file not found',0,.TRUE.)
     CALL write_line(trim(problem_volumes(volume_number)%filename),0,.TRUE.)
     STOP
  
END SUBROUTINE read_volume_list
