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
!3       ! number of excitations
!1       ! EXCITATION NUMBER
!impulse
!1.0 0.0
!2       ! EXCITATION NUMBER
!gaussian
!1.0 1e-9 3e-9
!3       ! EXCITATION NUMBER
!filename
!2   # number of header lines in file
!1   # column for time data
!2   # column for value data
!1.0 1e-9
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!     excitation function defined in a file: 23/9/2013
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
character*256	:: saved_input_line

integer	:: loop,sample,n_samples
real*8	:: t_in,ft_in
	   
integer   :: n_ignore
integer   :: n_time
integer   :: n_function
integer   :: n_max_read

real*8,allocatable	  :: read_data(:)

! START  

  CALL write_line('CALLED: read_excitation_function_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005,end=9005)n_excitation_functions
  
  CALL write_line_integer('number of excitation functions',n_excitation_functions,0,output_to_screen_flag)
  
  if ( allocated( excitation_functions ) ) GOTO 9000
  
  allocate ( excitation_functions(1:n_excitation_functions) )

  do excitation_number=1,n_excitation_functions
  
    CALL write_line_integer('Reading excitation number',excitation_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005,end=9005)read_number
    if (read_number.ne.excitation_number) goto 9010
 
! read excitation type string
    read(input_file_unit,'(A)',err=9010,end=9010)input_line
   
    CALL write_line( '...Reading excitation_function_type:'//trim(input_line),0,.TRUE. )
! save the original input line before converting to lower case - may be needed if it is a file name
    saved_input_line=input_line
    
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
   
    else if (input_line.eq.'differential_gaussian') then
    
      n_params=3
      excitation_functions(excitation_number)%type=excitation_function_type_differential_gaussian
   
    else if (input_line.eq.'noise') then
    
      n_params=1
      excitation_functions(excitation_number)%type=excitation_function_type_noise
   
    else if (input_line.eq.'sinusoidal_pulse') then
    
      n_params=7
      excitation_functions(excitation_number)%type=excitation_function_type_sinusoidal_pulse
   
    else if (input_line.eq.'ring_wave') then
    
      n_params=4
      excitation_functions(excitation_number)%type=excitation_function_type_ring_wave

    else

! we assume that the excitation type is a filename so try and open the file and read the contents

! go back to the original input line read from the input file and check whether this is a function in a file

      excitation_functions(excitation_number)%filename=saved_input_line

      write(*,*)'Excitation function file name='
      write(*,*)trim(excitation_functions(excitation_number)%filename)
      write(*,*)'Opening file:',trim(excitation_functions(excitation_number)%filename)

      open(UNIT=excitation_function_unit,						&
           FILE=excitation_functions(excitation_number)%filename,	&
           STATUS='old',							&
           ERR=9030)
	   
      write(*,*)'Enter the number of lines to ignore at the top of the data file'
      read(input_file_unit,*,err=9035,end=9035)n_ignore
      
      write(*,*)'Enter the column number for time data'
      read(input_file_unit,*,err=9035,end=9035)n_time
      
      write(*,*)'Enter the column number for function data'
      read(input_file_unit,*,err=9035,end=9035)n_function
      
      n_max_read=max(n_time,n_function)
      
      ALLOCATE( read_data(1:n_max_read) )
	   
      do loop=1,2

! read lines to ignore
        do i=1,n_ignore
          read(excitation_function_unit,*,end=1000,ERR=9040)
        end do
	
        sample=0
       
10      continue

          read(excitation_function_unit,*,end=1000,ERR=9040)(read_data(i),i=1,n_max_read)
	  sample=sample+1
	  if (loop.eq.2) then  
	    excitation_functions(excitation_number)%time_values_from_file(sample)=read_data(n_time)
	    excitation_functions(excitation_number)%function_values_from_file(sample)=read_data(n_function)
	  end if
	 
        goto 10  ! read next sample
        
1000    continue 

        n_samples=sample
        excitation_functions(excitation_number)%n_values_from_file=n_samples
	 
        if (loop.eq.1) then
	  allocate ( excitation_functions(excitation_number)%time_values_from_file(1:n_samples) )
	  allocate ( excitation_functions(excitation_number)%function_values_from_file(1:n_samples) )
	end if
	 
	rewind(unit=excitation_function_unit)
	 
      end do ! next loop
      
      excitation_functions(excitation_number)%type=excitation_function_type_file
                
      close(UNIT=excitation_function_unit)
      
      DEALLOCATE( read_data )
      
      write(*,*)'Number of excitation samples read:',excitation_functions(excitation_number)%n_values_from_file
      
      n_params=2

    end if
    
! read parameters
    excitation_functions(excitation_number)%n_parameters=n_params
    if (n_params.gt.0) then
      read(input_file_unit,*,err=9020,end=9020)(excitation_functions(excitation_number)%parameters(i),i=1,n_params)
    end if

  end do ! next excitation function

  CALL write_line('FINISHED: read_excitation_function_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating excitation_functions:',0,.TRUE.)
     CALL write_line('excitation_functions already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
  
9005 CALL write_line('Error reading excitation function list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9010 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Excitation functions should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9020 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Error reading parameters',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9030 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Unrecognised excitation function:',0,.TRUE.)
     write(*,*)trim(input_line)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9035 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Problem reading excitation function information:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
     
9040 CALL write_line('Error reading excitation function list packet',0,.TRUE.)
     CALL write_line('Problem reading excitation function from file:',0,.TRUE.)
     write(*,*)trim(input_line)
     CALL write_error_line(input_file_unit)
     STOP 1
     
  
END SUBROUTINE read_excitation_function_list
