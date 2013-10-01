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
! SUBROUTINE write_frequency_domain_data
! SUBROUTINE read_frequency_domain_data
!
! NAME
!    write_frequency_domain_data
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
SUBROUTINE write_frequency_domain_data(function_number)

USE post_process
USE file_information
USE file_header
USE output_formats

IMPLICIT NONE

  integer	:: function_number

! local variables

  character(len=256)	:: filename
  integer		:: frequency_loop

! START

  write(*,*)'Enter the frequency domain data output filename'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
    
  OPEN(unit=local_file_unit,file=filename)
  
  CALL write_frequency_domain_header_data(local_file_unit,1,function_of_frequency(function_number)%n_frequencies)

  do frequency_loop=1,function_of_frequency(function_number)%n_frequencies

    write(local_file_unit,frequency_domain_output_format)	&
       function_of_frequency(function_number)%frequency(frequency_loop),1,	&
       dble(function_of_frequency(function_number)%value(frequency_loop)),	&
       dimag(function_of_frequency(function_number)%value(frequency_loop)),	&
       function_of_frequency(function_number)%magnitude(frequency_loop),	&
       function_of_frequency(function_number)%phase(frequency_loop),	&
       function_of_frequency(function_number)%dB(frequency_loop)

  end do ! next frequency value
  
  CLOSE(unit=local_file_unit)

  RETURN
  
  
END SUBROUTINE write_frequency_domain_data
!
! NAME
!    read_frequency_domain_data
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
SUBROUTINE read_frequency_domain_data(function_number)

USE post_process
USE file_information
USE file_header
USE output_formats

IMPLICIT NONE

integer	:: function_number

! local variables

  character(len=256)	:: filename
  integer		:: n_frequencies
  integer		:: n_data_points
  integer		:: output_point
  
  integer		:: frequency
  integer		:: n_frequencies_read
  integer		:: loop
  
  real*8		:: f_in
  integer		:: n_in
  real*8		:: real_value_in
  real*8		:: imag_value_in
  real*8		:: mag_value_in
  real*8		:: phase_value_in
  real*8		:: dB_value_in

  character(len=256)	:: command
  
  logical		:: file_exists

! START

5  write(*,*)
  write(*,*)'Frequency domain output files:'
  
  command='ls -ltr *.fout'
  CALL system(command)

  write(*,*)'Enter the frequency domain filename'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
  
  CALL read_frequency_domain_header_data(local_file_unit,n_data_points,n_frequencies)
  
  write(*,*)'Opened file:',trim(filename)
  write(*,*)'Number of output points in file=',n_data_points
  write(*,*)'Number of frequencies in file=',n_frequencies
  
  write(*,*)'Enter the output point required'
  read(*,*)output_point
  write(record_user_inputs_unit,*)output_point

! read the input file   
       
  do loop=1,2
  
  frequency=0
  
  if (loop.EQ.2) then ! read header data again
    CALL read_frequency_domain_header_data(local_file_unit,n_data_points,n_frequencies)
  end if
  
10  CONTINUE

      read(local_file_unit,*,end=1000)f_in,n_in,real_value_in,imag_value_in,	&
                                      mag_value_in,phase_value_in,dB_value_in
      
      if (n_in.eq.output_point) then ! read this output
      
   	frequency=frequency+1
	
   	if (loop.eq.2) then
	
   	  function_of_frequency(function_number)%frequency(frequency)=f_in
   	  function_of_frequency(function_number)%value(frequency)=cmplx(real_value_in,imag_value_in)
   	  function_of_frequency(function_number)%magnitude(frequency)=mag_value_in
   	  function_of_frequency(function_number)%phase(frequency)=phase_value_in
   	  function_of_frequency(function_number)%dB(frequency)=dB_value_in
	  
   	end if ! loop.EQ.2
	
      end if ! n_in.eq.output_point so read this output
    
    GOTO 10  ! read next line of data
   
1000 CONTINUE

    n_frequencies_read=frequency
    
    if (loop.eq.1) then
      ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies_read) )
    end if
    
    rewind(unit=local_file_unit)
    
  end do ! next loop
  	 
  CLOSE(unit=local_file_unit)
  
  function_of_frequency(function_number)%n_frequencies=n_frequencies_read
  
  write(*,'(A,I10,A)')'Read ',n_frequencies_read,' frequency data values from file'

  RETURN
  
  
END SUBROUTINE read_frequency_domain_data
!
! NAME
!    read_frequency_domain_S_parameter_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     Maybe we should define a header for S parameter data...
!
! HISTORY
!
!     started 19/09/2013 CJS
!
!
SUBROUTINE read_frequency_domain_S_parameter_data(function_number1,function_number2)

USE post_process
USE file_information
USE file_header
USE output_formats

IMPLICIT NONE

integer	:: function_number1,function_number2

! local variables

  character(len=256)	:: filename
  integer		:: n_frequencies
  integer		:: n_data_points
  integer		:: output_point
  
  integer		:: frequency
  integer		:: n_frequencies_read
  integer		:: loop
  
  real*8		:: f_in
  
  real*8		:: real_value_S11_in
  real*8		:: imag_value_S11_in
  real*8		:: mag_value_S11_in
  real*8		:: dB_value_S11_in
  
  real*8		:: real_value_S21_in
  real*8		:: imag_value_S21_in
  real*8		:: mag_value_S21_in
  real*8		:: dB_value_S21_in

  character(len=256)	:: command
  
  logical		:: file_exists

! START

5  write(*,*)
  write(*,*)'Frequency domain output files:'
  
  command='ls -ltr *.fout'
  CALL system(command)

  write(*,*)'Enter the frequency domain filename'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
  
  read(local_file_unit,*) ! read header line
  
  write(*,*)'Opened file:',trim(filename)

! read the input file   
       
  do loop=1,2
  
    frequency=0
  
    if (loop.EQ.2) then ! read header line again
      read(local_file_unit,*) ! read header line
    end if
  
10  CONTINUE

      read(local_file_unit,*,end=1000)f_in,real_value_S11_in,imag_value_S11_in,mag_value_S11_in,dB_value_S11_in, &
                                           real_value_S21_in,imag_value_S21_in,mag_value_S21_in,dB_value_S21_in
                                          					         
      frequency=frequency+1
	
      if (loop.eq.2) then
	
   	function_of_frequency(function_number1)%frequency(frequency)=f_in
   	function_of_frequency(function_number1)%value(frequency)=cmplx(real_value_S11_in,imag_value_S11_in)
   	function_of_frequency(function_number1)%magnitude(frequency)=mag_value_S11_in
   	function_of_frequency(function_number1)%phase(frequency)=atan2(imag_value_S11_in,real_value_S11_in)
   	function_of_frequency(function_number1)%dB(frequency)=dB_value_S11_in
	
   	function_of_frequency(function_number2)%frequency(frequency)=f_in
   	function_of_frequency(function_number2)%value(frequency)=cmplx(real_value_S21_in,imag_value_S21_in)
   	function_of_frequency(function_number2)%magnitude(frequency)=mag_value_S21_in
   	function_of_frequency(function_number2)%phase(frequency)=atan2(imag_value_S21_in,real_value_S21_in)
   	function_of_frequency(function_number2)%dB(frequency)=dB_value_S21_in
	  
      end if ! loop.EQ.2
    
    GOTO 10  ! read next line of data
   
1000 CONTINUE

    n_frequencies_read=frequency
    
    if (loop.eq.1) then
      ALLOCATE ( function_of_frequency(function_number1)%frequency(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number1)%value(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number1)%magnitude(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number1)%phase(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number1)%dB(1:n_frequencies_read) )
      
      ALLOCATE ( function_of_frequency(function_number2)%frequency(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number2)%value(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number2)%magnitude(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number2)%phase(1:n_frequencies_read) )
      ALLOCATE ( function_of_frequency(function_number2)%dB(1:n_frequencies_read) )
    end if
    
    rewind(unit=local_file_unit)
    
  end do ! next loop
  	 
  CLOSE(unit=local_file_unit)
  
  function_of_frequency(function_number1)%n_frequencies=n_frequencies_read
  function_of_frequency(function_number2)%n_frequencies=n_frequencies_read
  
  write(*,'(A,I10,A)')'Read ',n_frequencies_read,' frequency data values from file'

  RETURN
  
  
END SUBROUTINE read_frequency_domain_S_parameter_data
