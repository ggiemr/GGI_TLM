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
! SUBROUTINE read_surface_list
!
! NAME
!     read_surface_list
!
! DESCRIPTION
!     read surface list packet
!
! Example packet:
!
!	Surface_list
!	3   Number of surfaces (integer)
!	1       surface_number (integer)
!	sphere 
!	1.0       surface parameters (n*real)
!	0.0 0.0 0.0
!	0.0 0.0 0.0
!	2       surface_number (integer)
!	rectangular_block
!	2.5 2.5 2.5      surface parameters (n*real)
!	0.0 0.0 0.0
!	0.0 0.0 0.0
!	3       surface_number (integer)
!	rectangular_block
!	3.0 3.0 3.0     surface parameters (n*real)
!	0.0 0.0 0.0
!	0.0 0.0 0.0
!
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
SUBROUTINE read_surface_list

USE TLM_general
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: surface_number
integer :: read_number

integer	:: n_params
integer :: i

character*256	:: input_line
character*256	:: tri_surface_filename

logical	:: file_exists

! START  

  CALL write_line('CALLED: Read_surface_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_surfaces
  
  CALL write_line_integer('number of surfaces',n_surfaces,0,output_to_screen_flag)
  
  if ( allocated( problem_surfaces ) ) GOTO 9000
  
  allocate ( problem_surfaces(1:n_surfaces) )

  do surface_number=1,n_surfaces
  
    CALL write_line_integer('Reading surface number',surface_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.surface_number) goto 9010
    
    problem_surfaces(surface_number)%surface_number=read_number
 
! read surface type string
    read(input_file_unit,'(A)',end=9010)input_line
   
    CALL write_line( '...STARTED reading surface_type:'//trim(input_line),0,.TRUE. )
    
! convert text to lower case
    CALL convert_to_lower_case(input_line,256)
   
    if (input_line.eq.'rectangular_block2') then
      n_params=6
      problem_surfaces(surface_number)%surface_type=surface_type_rectangular_block2

    else if (input_line.eq.'rectangular_block') then
      n_params=3
      problem_surfaces(surface_number)%surface_type=surface_type_rectangular_block
     
    else if (input_line.eq.'cylinder') then
      n_params=2
      problem_surfaces(surface_number)%surface_type=surface_type_cylinder

    else if (input_line.eq.'sphere') then
      n_params=1
      problem_surfaces(surface_number)%surface_type=surface_type_sphere
      
    else if (input_line.eq.'rectangle') then    
      n_params=2
      problem_surfaces(surface_number)%surface_type=surface_type_rectangle
     
    else if (input_line.eq.'circle') then
      n_params=1
      problem_surfaces(surface_number)%surface_type=surface_type_circle
     
    else if (input_line.eq.'quad') then
      n_params=12
      problem_surfaces(surface_number)%surface_type=surface_type_quad
     
    else if (input_line.eq.'xplane') then
      n_params=6
      problem_surfaces(surface_number)%surface_type=surface_type_xplane
     
    else if (input_line.eq.'yplane') then
      n_params=6
      problem_surfaces(surface_number)%surface_type=surface_type_yplane
     
    else if (input_line.eq.'zplane') then
      n_params=6
      problem_surfaces(surface_number)%surface_type=surface_type_zplane
     
    else if (input_line.eq.'triangle') then
      n_params=9
      problem_surfaces(surface_number)%surface_type=surface_type_triangle
     
    else if (input_line.eq.'triangulated_surface') then
    
      n_params=2
      problem_surfaces(surface_number)%surface_type=surface_type_triangulated_surface

! for a triangulated surface, read the mesh filename and check that the file exists
      
      read(input_file_unit,'(A)',end=9010)tri_surface_filename
      problem_surfaces(surface_number)%filename=tri_surface_filename
      
      CALL write_line('Checking the existance of file:',0,.TRUE.)
      CALL write_line(trim(problem_surfaces(surface_number)%filename),0,.TRUE.)
      
      inquire(file= problem_surfaces(surface_number)%filename,EXIST=file_exists)
      
      if (.NOT.file_exists) then
! error - no triangulated surface file exists
        goto 9050
      end if
     
    else if (input_line.eq.'vtk_triangulated_surface') then
    
      n_params=2
      problem_surfaces(surface_number)%surface_type=surface_type_vtk_triangulated_surface

! for a triangulated surface, read the mesh filename and check that the file exists
      
      read(input_file_unit,'(A)',end=9010)tri_surface_filename
      problem_surfaces(surface_number)%filename=tri_surface_filename
      
      CALL write_line('Checking the existance of file:',0,.TRUE.)
      CALL write_line(trim(problem_surfaces(surface_number)%filename),0,.TRUE.)
      
      inquire(file= problem_surfaces(surface_number)%filename,EXIST=file_exists)
      
      if (.NOT.file_exists) then
! error - no triangulated surface file exists
        goto 9060
      end if
	      
    else
    
      goto 9030

    end if
    
! read parameters
    problem_surfaces(surface_number)%n_surface_parameters=n_params
    if (n_params.gt.0) then
      read(input_file_unit,*,end=9020)    &
      (problem_surfaces(surface_number)%surface_parameters(i),i=1,n_params)
    end if
        
! read transformation
    problem_surfaces(surface_number)%trans%trans_type='euler'
    problem_surfaces(surface_number)%trans%trans_number=1
    read(input_file_unit,*,end=9040)    &
    (problem_surfaces(surface_number)%trans%parameters(i),i=1,6)
    
! convert euler angles to radians
    problem_surfaces(surface_number)%trans%parameters(1:3)=pi/180d0	&
       *problem_surfaces(surface_number)%trans%parameters(1:3)
   
    CALL write_line( '...FINISHED reading surface_type:'//trim(input_line),0,.TRUE. )
  
  end do ! next surface in surface_list


  CALL write_line('FINISHED: Read_surface_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating problem_surfaces:',0,.TRUE.)
     CALL write_line('problem_surfaces already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading surface_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading Surface_list_packet_data',0,.TRUE.)
     CALL write_line('surfaces should be numbered in order at the moment...',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading Surface_list_packet_data',0,.TRUE.)
     CALL write_line_integer('Number of paramters expected=',n_params,0,.TRUE.)
     STOP
  
9030 CALL write_line('Error reading Surface_list_packet_data',0,.TRUE.)
     CALL write_line('Unknown surface type:'//trim(input_line),0,.TRUE.)
     STOP
  
9040 CALL write_line('Error reading Surface_list_packet_data',0,.TRUE.)
     CALL write_line('Error reading transformation data',0,.TRUE.)
     STOP
  
9050 CALL write_line('Error in read_Surface_list_packet_data',0,.TRUE.)
     CALL write_line('Triangulated surface file not found',0,.TRUE.)
     CALL write_line(trim(problem_surfaces(surface_number)%filename),0,.TRUE.)
     STOP
  
9060 CALL write_line('Error in read_Surface_list_packet_data',0,.TRUE.)
     CALL write_line('Vtk format triangulated surface file not found',0,.TRUE.)
     CALL write_line(trim(problem_surfaces(surface_number)%filename),0,.TRUE.)
     STOP
  
END SUBROUTINE read_surface_list
