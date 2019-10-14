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
! SUBROUTINE write_time_domain_data_array
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
  
  integer               :: len_filename
  character(len=4)      :: extension

  integer		:: n_data_points
  integer		:: n_timesteps
  
  integer		:: output_point
  integer		:: loop
  integer		:: timestep
  
  real*8		:: t_in,value_in
  integer		:: n_in
  integer		:: n_timesteps_read

  integer               :: n_head,t_col,V_col,n_col

real*8,allocatable :: data_line(:)

  character(len=256)	:: command
  integer               :: status
  
  logical		:: file_exists
  logical               :: file_type_tout
  
  integer               :: i
  
! START

5 write(*,*)
  write(*,*)'Time domain output files:'
  
  command='ls -ltr *.tout'
  CALL system(command,status)
  
  if (status.NE.0) then       ! no tout files were found so list all files
    write(*,*)
    write(*,*)' No .tout files found, listing all files...'
    write(*,*)    
    command='ls -ltr '
    CALL system(command,status)
  end if

  write(*,*)'Enter the time domain filename'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist:',trim(filename)
    STOP
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  write(post_process_info_unit,*)'	Time domain data filename:',trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
  write(*,*)'Opened file:',trim(filename)
  
! Check whether this is a file with extension tout. If it is, read in the normal way, otherwise read a more general format

  len_filename=len(trim(filename))
  if (len_filename.GE.4) then
    extension=filename(len_filename-3:len_filename)
  else
    extension='****'   ! something that is not 'tout'
  end if
  
  if (extension.EQ.'tout') then
  
    CALL read_time_domain_header_data(local_file_unit,n_data_points,n_timesteps)
  
    write(*,*)'Number of output points in file=',n_data_points
    write(*,*)'Number of timesteps in file=',n_timesteps
  
    write(*,*)'Enter the output point required'
    read(*,*)output_point
    write(record_user_inputs_unit,*)output_point
  
    write(post_process_info_unit,*)'	Time domain output point:',output_point
        
    rewind(unit=local_file_unit)
    
    file_type_tout=.TRUE.
    
    n_head=3
    t_col=1
    V_col=3
    
  else
   
    write(*,*)'Enter the number of header lines in the file'
    read(*,*)n_head
    write(record_user_inputs_unit,*)n_head,'  # number of header lines in the file'

    write(*,*)'Enter the column number for time data'
    read(*,*)t_col
    write(record_user_inputs_unit,*)t_col,'  # time data column'

    write(*,*)'Enter the column number for voltage data'
    read(*,*)V_col
    write(record_user_inputs_unit,*)V_col,'  # Voltage data column'
    
    file_type_tout=.FALSE.
    output_point=1
    
  end if
  
! read the input file   

  n_col=max(t_col,V_col)

  ALLOCATE( data_line(1:n_col) )
       
  do loop=1,2
  
! read file header
    do i=1,n_head
      read(local_file_unit,*)
    end do
  
    timestep=0
    
10  CONTINUE

      if (file_type_tout) then
      
        read(local_file_unit,*,end=1000)t_in,n_in,value_in  ! tout format
        
      else 
      
        read(local_file_unit,*,END=1000)(data_line(i),i=1,n_col)
        t_in =data_line(t_col)
        value_in=data_line(V_col)
        n_in=1
        
      end if
      
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
  
  DEALLOCATE( data_line )
  
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
    
  write(post_process_info_unit,*)'	Time domain output data filename:',trim(filename)

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
!
! NAME
!    write_time_domain_data_array
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/11/2013 CJS
!
!
SUBROUTINE write_time_domain_data_array(n_functions)

USE post_process
USE file_information
USE file_header
USE output_formats

IMPLICIT NONE

integer			:: n_functions

! local variables

integer			:: function_number

  character(len=256)	:: filename

  integer		:: n_data_points
  integer		:: n_timesteps
  
  integer		:: timestep
    
! START
  
  write(*,*)'Enter the filename for the time domain data'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
    
  write(post_process_info_unit,*)'	Time domain output array data filename:',trim(filename)

  OPEN(unit=local_file_unit,file=filename)
  
  n_timesteps  =function_of_time(1)%n_timesteps
  
  n_data_points=n_functions
  
  CALL write_time_domain_header_data(local_file_unit,n_data_points,n_timesteps)
    
! write the data file   

  do function_number=1,n_functions
       
    n_timesteps  =function_of_time(function_number)%n_timesteps
    
    do timestep=1,n_timesteps
  
      write(local_file_unit,time_domain_output_format)	&
                                      function_of_time(function_number)%time(timestep),function_number,	&
                                      function_of_time(function_number)%value(timestep)
          
    end do ! next timestep
    
  end do ! next function
  	 
  CLOSE(unit=local_file_unit)

  RETURN
  
  
END SUBROUTINE write_time_domain_data_array
