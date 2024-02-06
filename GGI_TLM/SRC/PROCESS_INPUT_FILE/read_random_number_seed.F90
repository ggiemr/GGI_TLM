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
! SUBROUTINE read_random_number_seed
!
! NAME
!     read_random_number_seed
!
! DESCRIPTION
!     read a seed value for the random number generator and 
!     reset the seed value
!     Randomising the seed value is based on the subroutine on the webpage:
!
!     http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!
! Example packet:
!
!
!random_number_seed
!3   ! number of values
!12349876
!54723856
!78364562
!
! COMMENTS
!     
!
! HISTORY
!
!     started 29/01/2014 CJS
!
!
SUBROUTINE read_random_number_seed

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

  integer :: nvalues
  real*8  :: r_random
  
  integer, allocatable :: seed(:)
  integer	:: i,n
  
! START  

  CALL write_line('CALLED: read_random_number_seed',0,output_to_screen_flag)

! read reflection coefficient data

  call random_seed(size = n)
  
  write(*,*)'Number of random number seed values required=',n
  
  read(input_file_unit,*,err=9000,end=9000)nvalues
    
  if (nvalues.eq.0) then
  
    CALL set_random_seed()
    
  else if (nvalues.lt.n) then
! Error: insufficient number of values for for the random_seed function  
   GOTO 9010
  
  else
    
    ALLOCATE(seed(n))
    
! read the first n values from the input file  
    do i=1,n   
      read(input_file_unit,*,err=9000)seed(i)
    end do
  
    call random_seed(put=seed)

    DEALLOCATE(seed)

! evaluate a few random numbers to 'warm up' the random number generator if small values of seed are used for example  
    do i=1,256
      CALL random_number(r_random)
    end do

  end if  ! nvalues.eq.0
  
  CALL write_line('FINISHED: read_random_number_seed',0,output_to_screen_flag)
       
  RETURN
  
9000 CALL write_line('Error reading random_number_seed packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
  
9010 CALL write_line('Error reading random_number_seed packet data from input file:',0,.TRUE.)
     CALL write_line_integer('Number of values should be at least ',n,0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
       
END SUBROUTINE read_random_number_seed
