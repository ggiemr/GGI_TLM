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
! SUBROUTINE read_time_domain_data
! SUBROUTINE write_time_domain_data
!
! NAME
!    read_time_domain_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE read_time_domain_data(function_number)

USE post_process
USE file_information

IMPLICIT NONE

integer	:: function_number

! local variables

character(len=256)	:: filename

  integer		:: n_data_points
  integer		:: n_timesteps
  
  integer		:: output_point
  integer		:: loop
  integer		:: timestep
  
  real*8		:: t_in,value_in
  integer		:: n_in
  integer		:: n_timesteps_read

  character(len=256)	:: command
  
  logical		:: file_exists
  
! START

5 write(*,*)
  write(*,*)'Time domain output files:'
  
  command='ls -ltr *.tout'
  CALL system(command)

  write(*,*)'Enter the time domain filename'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
  
  CALL read_time_domain_header_data(local_file_unit,n_data_points,n_timesteps)
  
  write(*,*)'Opened file:',trim(filename)
  write(*,*)'Number of output points in file=',n_data_points
  write(*,*)'Number of timesteps in file=',n_timesteps
  
  write(*,*)'Enter the output point required'
  read(*,*)output_point
  write(record_user_inputs_unit,*)output_point
  
! read the input file   
       
  do loop=1,2
  
  timestep=0
  
  if (loop.EQ.2) then ! read header data again
    CALL read_time_domain_header_data(local_file_unit,n_data_points,n_timesteps)
  end if
  
10  CONTINUE

      read(local_file_unit,*,end=1000)t_in,n_in,value_in
      
      if (n_in.eq.output_point) then ! read this output
      
   	timestep=timestep+1
	
   	if (loop.eq.2) then
	
   	  function_of_time(function_number)%time(timestep)=t_in
   	  function_of_time(function_number)%value(timestep)=value_in
	  
   	end if ! loop.EQ.2
	
      end if ! n_in.eq.output_point so read this output
    
    GOTO 10  ! read next line of data
   
1000 CONTINUE

    n_timesteps_read=timestep
    
    if (loop.eq.1) then
      ALLOCATE ( function_of_time(function_number)%time(1:n_timesteps_read) )
      ALLOCATE ( function_of_time(function_number)%value(1:n_timesteps_read) )
    end if
    
    rewind(unit=local_file_unit)
    
  end do ! next loop
  	 
  CLOSE(unit=local_file_unit)
  
  function_of_time(function_number)%n_timesteps=n_timesteps_read
  
  write(*,'(A,I10,A)')'Read ',n_timesteps_read,' timestep data values from file'

  RETURN
  
  
END SUBROUTINE read_time_domain_data
!
! NAME
!    write_time_domain_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE write_time_domain_data(function_number)

USE post_process
USE file_information
USE file_header
USE output_formats

IMPLICIT NONE

integer			:: function_number

! local variables

  character(len=256)	:: filename

  integer		:: n_data_points
  integer		:: n_timesteps
  
  integer		:: timestep
    
! START
  
  write(*,*)'Enter the filename for the time domain curve'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)

  OPEN(unit=local_file_unit,file=filename)
  
  n_data_points=1
  n_timesteps  =function_of_time(function_number)%n_timesteps
  
  CALL write_time_domain_header_data(local_file_unit,n_data_points,n_timesteps)
    
! write the input file   
       
  do timestep=1,n_timesteps
  
    write(local_file_unit,time_domain_output_format)	&
                                     function_of_time(function_number)%time(timestep),1,	&
                                     function_of_time(function_number)%value(timestep)
          
  end do ! next timestep
  	 
  CLOSE(unit=local_file_unit)

  RETURN
  
  
END SUBROUTINE write_time_domain_data
