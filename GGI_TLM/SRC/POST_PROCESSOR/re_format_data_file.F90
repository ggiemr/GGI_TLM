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
! SUBROUTINE re_format_data_file
!
! NAME
!    re_format_data_file
!
! DESCRIPTION
!    Processes for re-formatting data files so that csv files from measurement can be converted into files 
!    suitable for plottting with gnuplot for example.
!   
!    The following processes are implemented here:
!    1. Strip file header data
!    2. Read csv format data
!    3. Re-arrange columns of data
!    4. Multiply a column by a given factor (eg scale from Hz to GHz)
!    5. Keep data within a given parameter range
!    6. Combine files
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2014 CJS
!     
!
SUBROUTINE re_format_data_file(function_number)

USE post_process
USE file_information

IMPLICIT NONE

integer	:: function_number

! local variables

  character(len=256)	:: filename

  integer		:: sample,n_samples
  
  logical		:: file_exists
	   
  integer   		:: first_line
  integer   		:: last_line
  
  integer   		:: n_lines
  
  integer   		:: n_columns
  
  integer   		:: n_add
  integer   		:: n_rescale
  
  integer,allocatable	:: column_list(:)
  integer		:: max_column

  real*8,allocatable	:: input_line_data(:)
  real*8,allocatable	:: data(:,:)
  real*8,allocatable	:: max_data(:)
  real*8,allocatable	:: min_data(:)
  real*8,allocatable	:: data_scale_factor(:)
  real*8,allocatable	:: data_add_constant(:)
  
  logical	:: output_range_check
  integer   :: range_check_column
  real*8    :: max_output
  real*8    :: min_output
  
  integer 		:: loop,i,col_loop
  
  character		:: ch
  
  character(len=256)	:: command

! START
  
  write(*,*)
  write(*,*)'Output files:'
  command='ls -ltr '
  CALL system(command)

! get the filename for the data
5 CONTINUE
  write(*,*)'Enter the filename for the input data'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
! open and read the file
  
  OPEN(unit=local_file_unit,file=filename)
  
  CALL write_file_format_information(local_file_unit,n_lines)
  
  rewind(unit=local_file_unit)
  
  write(*,*)'Enter the first line of the data file to process'
  read(*,*)first_line
  write(record_user_inputs_unit,*)first_line,' First line of the data file to process'
  
  write(*,*)'Enter the first line of the data file to process or 0 to read the whole file'
  read(*,*)last_line
  write(record_user_inputs_unit,*)last_line,' Last line of the data file to process or 0 to read the whole file'
  
  if ( (last_line.EQ.0).OR.(last_line.GT.n_lines) ) then
    last_line=n_lines
  end if 
  
  if (last_line.LT.first_line) then
    write(*,*)'Error: last line specified is less than the first line'
    STOP
  end if 
    
  write(*,*)'Enter the number of columns of data to read'
  read(*,*)n_columns
  write(record_user_inputs_unit,*)n_columns,' Number of columns of data to read'
    
  max_column=0
  ALLOCATE( column_list(1:n_columns) )
  
  do loop=1,n_columns
  
    write(*,'(A,I2,A)')'Which Column should output file column ',loop,' come from? '
    read(*,*)column_list(loop)
    write(record_user_inputs_unit,*)column_list(loop),' Column number for output file column',loop
    max_column=max(max_column,column_list(loop))
    
  end do 
  
! Initial read of the data file
    
  ALLOCATE( input_line_data(1:max_column) )
       
  do loop=1,2

! read lines to ignore
    do i=1,first_line-1
      read(local_file_unit,*,ERR=9000)
    end do
    
    n_samples=last_line-first_line+1
    
    do sample=1,n_samples
    
      read(local_file_unit,*,ERR=9000)(input_line_data(i),i=1,max_column)
      
      if (loop.eq.2) then  
! arrange the data into the correct columns for the output file

        do col_loop=1,n_columns
	
          data(sample,col_loop)=input_line_data( column_list(col_loop) )
	  
	  max_data(col_loop)=max( max_data(col_loop),data(sample,col_loop) )
	  min_data(col_loop)=min( min_data(col_loop),data(sample,col_loop) )
	  
	end do
	
      end if
     
    end do ! next sample to read
     
    if (loop.eq.1) then
    
      ALLOCATE ( data(1:n_samples,1:n_columns) )
      ALLOCATE ( max_data(1:n_columns) )
      ALLOCATE ( min_data(1:n_columns) )
      
      max_data(1:n_columns)=-1D30
      min_data(1:n_columns)= 1D30
      
    end if
     
    rewind(unit=local_file_unit)
     
  end do ! next loop
      
  close(UNIT=local_file_unit)

! Write column data ranges to the screen      
  write(*,*)'Number of data samples read:',n_samples
  write(*,*)'Data column ranges:'
  do col_loop=1,n_columns
    write(*,'(A,I3,A,E16.6,A,E16.6)')'Column ',col_loop,' min value=',min_data(col_loop),' max value=',max_data(col_loop)
  end do
  
! Restrict the data ranges if required     
  
! ADD CONSTANT TO DATA  
  ALLOCATE ( data_add_constant(n_columns) )
  data_add_constant(:)=0d0
    
  write(*,*)'Enter the number of columns of data to add constant to'
  read(*,*)n_add
  write(record_user_inputs_unit,*)n_add,' Number of columns of data to add constant to'
  
  do loop=1,n_add
  
    write(*,*)'Enter the column to add constant to'
    read(*,*)i
    write(record_user_inputs_unit,*)i,' Number of column to add constant to'
  
    write(*,*)'Constant to add'
    read(*,*)data_add_constant(i)
    write(record_user_inputs_unit,*)data_add_constant(i),' Constant to add'
  
  end do

! Apply constant addition 
  do col_loop=1,n_columns
    data(1:n_samples,col_loop)=data(1:n_samples,col_loop)+data_add_constant(col_loop)
  end do
  min_data(:)=min_data(:)+data_add_constant(:)
  max_data(:)=max_data(:)+data_add_constant(:)
  
  write(*,*)'Data column ranges:'
  do col_loop=1,n_columns
    write(*,'(A,I3,A,E16.6,A,E16.6)')'Column ',col_loop,' min value=',min_data(col_loop),' max value=',max_data(col_loop)
  end do
  
! RE_SCALE DATA  
  ALLOCATE ( data_scale_factor(n_columns) )
  data_scale_factor(:)=1d0
    
  write(*,*)'Enter the number of columns of data to re-scale'
  read(*,*)n_rescale
  write(record_user_inputs_unit,*)n_rescale,' Number of columns of data to re-scale'
  
  do loop=1,n_rescale
  
    write(*,*)'Enter the column to re-scale'
    read(*,*)i
    write(record_user_inputs_unit,*)i,' Number of column to re-scale'
  
    write(*,*)'Re-scale factor'
    read(*,*)data_scale_factor(i)
    write(record_user_inputs_unit,*)data_scale_factor(i),' Re-scale factor'
  
  end do

! Apply re-scaling factors  
  do col_loop=1,n_columns
    data(1:n_samples,col_loop)=data(1:n_samples,col_loop)*data_scale_factor(col_loop)
  end do
  min_data(:)=min_data(:)*data_scale_factor(:)
  max_data(:)=max_data(:)*data_scale_factor(:)
  
  write(*,*)'Data column ranges:'
  do col_loop=1,n_columns
    write(*,'(A,I3,A,E16.6,A,E16.6)')'Column ',col_loop,' min value=',min_data(col_loop),' max value=',max_data(col_loop)
  end do
  
  write(*,*)'Do you want to output data over a restricted range ? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A,A)')ch,' ! Restrict range of data?'
  
  if ( (ch.EQ.'y').OR.(ch.EQ.'Y') ) then

    output_range_check=.TRUE.
  
    write(*,*)'Enter the column for the range check'
    read(*,*)range_check_column
    write(record_user_inputs_unit,*)range_check_column,' Column for the range check'
  
    write(*,*)'Enter the minimum value of the output range'
    read(*,*)min_output
    write(record_user_inputs_unit,*)min_output,' Minimum value of the output range'
  
    write(*,*)'Enter the maximum value of the output range'
    read(*,*)max_output
    write(record_user_inputs_unit,*)max_output,' Maximum value of the output range'
      
  else

    output_range_check=.FALSE.
    
  end if

! open a file for the re-formatted data
  
  write(*,*)'Enter the filename for the formatted data'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
 
! Give the option to append to existing file...  
!! open the file to write
!  
!  OPEN(unit=local_file_unit,file=filename)
    
  write(*,*)'Do you want to append new data to an existing data file? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A,A)')ch,' ! Append to existing file?'
  
  if ( (ch.EQ.'y').OR.(ch.EQ.'Y') ) then

! open file and read to end of file if it exists
    inquire(file=trim(filename),exist=file_exists)
    if (file_exists) then
    
      OPEN(unit=local_file_unit,file=filename,status='OLD',access='APPEND',position='APPEND')
      
    else
    
! file does not exist so open new file    
      OPEN(unit=local_file_unit,file=filename)
 
    end if

  else  

! open new file 
    OPEN(unit=local_file_unit,file=filename)
    
  end if

  do sample=1,n_samples       
    
! Eliminate absolute data values smaller than 1E-30 - smaller numbers cause format problems
    do col_loop=1,n_columns
      if ( abs(data(sample,col_loop) ).LT.1D-30) data(sample,col_loop)=0d0
    end do

    if (output_range_check) then
    
! do range test
      if (  (data(sample,range_check_column).GE.min_output).AND.(data(sample,range_check_column).LE.max_output) ) then
 ! Write column data to file    
        write(local_file_unit,*)(data(sample,col_loop),col_loop=1,n_columns)
      end if
      
    else

! Write column data to file    
      write(local_file_unit,*)(data(sample,col_loop),col_loop=1,n_columns)
      
    end if
    
  end do
  
  close(unit=local_file_unit)
  
  DEALLOCATE( column_list )
  DEALLOCATE( input_line_data )
  DEALLOCATE( data )
  DEALLOCATE( min_data )
  DEALLOCATE( max_data ) 
  DEALLOCATE( data_scale_factor )
  DEALLOCATE( data_add_constant )

  RETURN
     
9000 CALL write_line('Error reading data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(filename)
     STOP 
  
END SUBROUTINE re_format_data_file
!
! _________________________________________________________________________________
!
!
SUBROUTINE write_file_format_information(unit,number_of_lines)

IMPLICIT NONE

integer :: unit
integer :: number_of_lines

! local variables

character*256 input_line

integer :: number_of_data_lines
integer :: number_of_comment_lines

logical :: was_last_line_a_number
logical :: is_line_a_number

integer :: block_start_line
integer :: block_end_line
integer :: block_length

! function variables

logical :: is_a_number

! START

write(*,*)' '
write(*,*)'Summary of data file structure'
write(*,*)' '

number_of_lines=0
number_of_data_lines=0
number_of_comment_lines=0

! read first line and work out if this is a number (data line) or not

  read(unit,'(A)',END=1000)input_line
  number_of_lines=number_of_lines+1
  was_last_line_a_number=is_a_number(input_line)
  block_start_line=1

! add this line to the appropriate count    
  if (was_last_line_a_number) then
    number_of_data_lines=number_of_data_lines+1
  else
    number_of_comment_lines=number_of_comment_lines+1  
  end if

10 CONTINUE
    read(unit,'(A)',END=1000)input_line
    
    number_of_lines=number_of_lines+1

! check whether this line is a number or comment line    
    is_line_a_number=is_a_number(input_line)

! add this line to the appropriate count    
    if (is_line_a_number) then
      number_of_data_lines=number_of_data_lines+1
    else
      number_of_comment_lines=number_of_comment_lines+1  
    end if

! see if we have come to the end of a block of numbers or comments 
   
    if (is_line_a_number.NEQV.was_last_line_a_number) then
! the current line type is different to the previous line type      
    
      block_end_line=number_of_lines-1
      block_length=block_end_line-block_start_line+1
            
      if (was_last_line_a_number) then
! Numeric data block    
        if (block_length.EQ.1) then
          write(*,*)block_length,' Line of numeric data, line number',block_start_line
        else
          write(*,*)block_length,' Lines of numeric data, start line',block_start_line,' end line',block_end_line
        end if
	
      else
! Comment block         
        if (block_length.EQ.1) then
          write(*,*)block_length,' Line of comment data, line number',block_start_line
        else
          write(*,*)block_length,' Lines of comment data, start line',block_start_line,' end line',block_end_line
        end if
      
      end if
! start of new block is the current line      
      block_start_line=number_of_lines
    
    end if
    
    was_last_line_a_number=is_line_a_number

    GOTO 10

1000 CONTINUE

! Final block
  block_end_line=number_of_lines
  block_length=block_end_line-block_start_line+1
  	
  if (was_last_line_a_number) then
! Numeric data block    
    if (block_length.EQ.1) then
      write(*,*)block_length,' Line of numeric data, line number',block_start_line
    else
      write(*,*)block_length,' Lines of numeric data, start line',block_start_line,' end line',block_end_line
    end if

  else
! Comment block         
    if (block_length.EQ.1) then
      write(*,*)block_length,' Line of comment data, line number',block_start_line
    else
      write(*,*)block_length,' Lines of comment data, start line',block_start_line,' end line',block_end_line
    end if
  
  end if

  write(*,*)'Number of lines        =',number_of_lines
  write(*,*)'Number of data lines   =',number_of_data_lines
  write(*,*)'Number of comment lines=',number_of_comment_lines

  write(*,*)' '

  RETURN

END SUBROUTINE write_file_format_information
!
! _________________________________________________________________________________
!
!
logical FUNCTION is_a_number(string)

IMPLICIT NONE

character*256 string

! local variables

  real*8	:: number
  integer	  :: stat

! START

  is_a_number=.FALSE.
  
  read(string,*,iostat=stat)  number
  
  if (stat.EQ.0) then
    is_a_number=.TRUE.
  else
    is_a_number=.FALSE.
  end if
  
!  write(*,*)'TEST_STRING:',trim(string),' TEST:',is_a_number

END  FUNCTION is_a_number
