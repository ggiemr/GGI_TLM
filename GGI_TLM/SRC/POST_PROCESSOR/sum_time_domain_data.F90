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
! SUBROUTINE sum_time_domain_data
!
! NAME
!    sum_time_domain_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/11/2014 CJS
!
!
SUBROUTINE sum_time_domain_data


USE post_process
USE file_information

IMPLICIT NONE

! local variables
  
  real*8	:: t1,t2,result
  
  integer	:: function_number
  integer	:: n_functions
  integer	:: n_timesteps
  integer	:: time_loop
  
  real*8,allocatable	:: multiplication_factor(:)
  
! START

!  write(*,*)'Sum Time Domain Data'
  
  write(*,*)' '
  write(*,*)'Result=sum 1..N {fn*multiplication_factor}'
  write(*,*)' '
  write(*,*)'where fn is a set of N time domain quantities'
  write(*,*)' '
  
  write(*,*)'Enter the number of time domain quantities to sum:'
  read(*,*)n_functions
  write(record_user_inputs_unit,*)n_functions,' number of functions to sum'
  
  ALLOCATE( multiplication_factor(1:n_functions) )
  
  write(post_process_info_unit,*)'	Number of functions to sum=',n_functions

  n_functions_of_time=n_functions+1
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  write(post_process_info_unit,*)'	Time domain functions:'
  do function_number=1,n_functions
  
    write(*,*)'File for time domain data:'
    CALL read_time_domain_data(function_number)
    
    write(*,*)'Enter the multiplication_factor for this function'
    read(*,*)multiplication_factor(function_number)
    write(record_user_inputs_unit,*)multiplication_factor(function_number),' Multiplication_factor for this function'

   write(post_process_info_unit,*)'	Multiplication_factor for this function:',multiplication_factor(function_number)
 
    if (function_number.gt.1) then
! check that the time samples match...
      if (function_of_time(1)%n_timesteps.NE.function_of_time(function_number)%n_timesteps) then
        write(*,*)'Time serioes mismatch in functions f 1 and f',function_number
        write(*,*)'n_timesteps, f1=',function_of_time(1)%n_timesteps
        write(*,*)'n_timesteps, fn=',function_of_time(function_number)%n_timesteps
        STOP
      end if
  
      do time_loop=1,function_of_time(1)%n_timesteps
  
        if ( function_of_time(1)%time(time_loop).NE.	&
             function_of_time(function_number)%time(time_loop) ) then
	 
          write(*,*)'time mismatch in functions f 1 and f',function_number
          write(*,*)'time number',time_loop
          write(*,*)'time, f1=',function_of_time(1)%time(time_loop)
          write(*,*)'time, fn=',function_of_time(function_number)%time(time_loop)
          STOP
      
        end if
  
      end do
      
    end if !function_number.gt.1
    
  end do ! next function number

! Allocate memory for result

  n_timesteps=function_of_time(1)%n_timesteps
  function_number=n_functions_of_time
  
  function_of_time(function_number)%n_timesteps=n_timesteps
  
  ALLOCATE ( function_of_time(function_number)%time(1:n_timesteps) )
  function_of_time(function_number)%time(1:n_timesteps)=		&
                function_of_time(1)%time(1:n_timesteps)
     
  ALLOCATE ( function_of_time(function_number)%value(1:n_timesteps) )
  
  do time_loop=1,n_timesteps
  
    result=(0d0,0d0)
    
    do function_number=1,n_functions

      result=result+function_of_time(function_number)%value(time_loop)*multiplication_factor(function_number)
      
    end do
    
    function_number=n_functions_of_time
    
    function_of_time(function_number)%value(time_loop)=result
    
  end do ! next time value
  
  function_number=n_functions_of_time
  CALL write_time_Domain_Data(function_number)
  
  DEALLOCATE( multiplication_factor )

  CALL Deallocate_post_data()

  RETURN
  

  
END SUBROUTINE sum_time_domain_data
