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
! SUBROUTINE read_data
!
! NAME
!     read data
!
! DESCRIPTION
!     read the complex frequency domain data 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 12/12/2012 CJS
!
!
SUBROUTINE read_data

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE
 
! local variables

  integer	:: read_loop
  integer	:: function_loop
  integer	:: freq_loop
  
  integer	:: n_samples(1:4)
  character(len=256) :: filename(1:4)
  
  real*8,allocatable		:: local_frequency(:,:)

  real*8	:: f_in,real_value_in,imag_value_in

! START

  CALL write_line('CALLED: read_data',0,ff_output_to_screen)

  if (fit_type.eq.dielectric_material) then 
  
    n_functions=1
    filename(1)=trim(FF_name)//dielectric_initial_data_extension
 
  else if (fit_type.eq.magnetic_material) then 
  
    n_functions=1
    filename(1)=trim(FF_name)//magnetic_initial_data_extension
 
  else if (fit_type.eq.thin_layer) then 
  
    n_functions=4
    filename(1)=trim(FF_name)//thin_layer_z11_initial_data_extension
    filename(2)=trim(FF_name)//thin_layer_z12_initial_data_extension
    filename(3)=trim(FF_name)//thin_layer_z21_initial_data_extension
    filename(4)=trim(FF_name)//thin_layer_z22_initial_data_extension
 
  else if (fit_type.eq.impedance) then 
   
    n_functions=1
    filename(1)=trim(FF_name)//impedance_initial_data_extension

  end if

! READ FILES IN TWO STAGES, FIRST GET THE NUMBER OF SAMPLES THEN ALLOCATE ARRAYS AND READ THE DATA
 
  do read_loop=1,2
  
    do function_loop=1,n_functions
  
      OPEN(unit=initial_data_file_unit,file=filename(function_loop),status='OLD',err=9000)
      n_samples(function_loop)=0
      
10    CONTINUE

      read(initial_data_file_unit,*,end=1000)f_in,real_value_in,imag_value_in
      
      n_samples(function_loop)=n_samples(function_loop)+1

      if (read_loop.eq.2) then

        local_frequency(function_loop,n_samples(function_loop))=f_in
        value	 (function_loop,n_samples(function_loop))=dcmplx(real_value_in,imag_value_in)
        
      end if ! read_loop.EQ.2
    
      GOTO 10  ! read next line of data
   
1000 CONTINUE

    CLOSE(unit=initial_data_file_unit)
    
  end do! next function

! check that the number of samples is the same for all files  
  
    n_values=n_samples(1)
    do function_loop=2,n_functions
      if(n_samples(function_loop).ne.n_values) GOTO 9010
    end do ! next function
    
    if (read_loop.eq.1) then
    
      ALLOCATE (  local_frequency(1:n_functions,n_values )  )
      ALLOCATE (  value    (1:n_functions,n_values )  )
      
    end if
    
  end do ! next read_loop

! create a single frequency list for all functions
  ALLOCATE (  frequency(1:n_values)  )
  
  frequency(1:n_values)=local_frequency(1,1:n_values)
  
! Check that the frequencies are the same in all data sets  
  do freq_loop=1,n_values
    do function_loop=1,n_functions
      if(local_frequency(function_loop,freq_loop).ne.frequency(freq_loop)) GOTO 9020
    end do
  end do ! next function

! Set the frequency normalisation constant to be the last   
  
  fnorm=frequency(n_values)
  wnorm=2d0*pi*fnorm

! construct the normalised frequency array  
  ALLOCATE (  normalised_frequency(1:n_values)  )
  normalised_frequency(1:n_values)=frequency(1:n_values)/fnorm
  
! construct the normalised complex frequency array  
  ALLOCATE (  s(1:n_values)  )
  s(1:n_values)=2d0*pi*j*frequency(1:n_values)/wnorm
    
  if (allocated( local_frequency )) DEALLOCATE ( local_frequency )
  
  CALL write_line_integer('Number of values read=',n_values,0,ff_output_to_screen)

  CALL write_line('FINISHED: read_data',0,ff_output_to_screen)
  
  RETURN
  
9000 CALL write_line('Error in read_data',0,.TRUE.)
     CALL write_line('Unable to open file:'//trim(filename(function_loop)),0,.TRUE.)
     STOP
  
9010 CALL write_line('Error in read_data',0,.TRUE.)
     CALL write_line('Number of values is not the same across all the files',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error in read_data',0,.TRUE.)
     CALL write_line('Maximum frequency is not the same across all the files',0,.TRUE.)
     STOP
       
END SUBROUTINE read_data

