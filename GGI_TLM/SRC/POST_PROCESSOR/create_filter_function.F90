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
! SUBROUTINE create_filter_function
!
! NAME
!    create_filter_function
!
! DESCRIPTION
!    create a low pass or high pass filter function
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 27/2/2018 CJS
!
!
SUBROUTINE create_filter_function

USE constants
USE post_process
USE filter_types		        
USE filter_functions		        
USE file_information

IMPLICIT NONE

! local variables

  integer			:: function_number

! Filter stuff
  character*256			:: filter_file_name
  type(Sfilter) 		:: Sfilter0
  type(Sfilter) 		:: Sfilter1
  type(Sfilter) 		:: Sfilter2
  type(Sfilter) 		:: Sfilter3
  type(Sfilter) 		:: Sfilter4

  real*8			:: frequency
  integer                       :: order

  character                     :: ch

! START

  write(*,*)'Please enter the type of filter required (LPF, HPF)'
  read(*,'(A)')ch

  write(*,*)'Please enter the cutoff frequency (Hz)'
  read(*,*)frequency

  write(*,*)'Please enter the order of the filter'
  read(*,*)order
  
  if (order.LT.1) then
    write(*,*)'The filter order should be greater than zero'
    STOP
  end if
  if (order.GT.10) then
    write(*,*)'The filter order should be less than or equal to 10'
    STOP
   end if
  
  if ( (ch.EQ.'L').OR.(ch.EQ.'l') ) then
! Low pass filter
    
    Sfilter0=Allocate_Sfilter(0,order)
    
    Sfilter0%wnorm=2d0*pi*frequency
    
    Sfilter0%a%coeff(0)=1d0

  else if ( (ch.EQ.'H').OR.(ch.EQ.'h') ) then
! High pass filter
    
    Sfilter0=Allocate_Sfilter(order,order)
    
    Sfilter0%wnorm=2d0*pi*frequency
    
    Sfilter0%a%coeff(0:order-1)=0d0
    Sfilter0%a%coeff(order)=1d0

  else
    write(*,*)'The filter type should be LPF or HPF'
    STOP
  end if
  
  Sfilter0%b%coeff(0)=1d0
  Sfilter0%b%coeff(order)=1d0
  
  if (order.EQ.2) then

    Sfilter0%b%coeff(1)=1.414d0 
 
  else if (order.EQ.3) then

    Sfilter0%b%coeff(1)=2d0 
    Sfilter0%b%coeff(2)=2d0 
 
  else if (order.EQ.4) then

    Sfilter0%b%coeff(1)=2.613d0 
    Sfilter0%b%coeff(2)=3.414d0 
    Sfilter0%b%coeff(3)=2.613d0 
 
  else if (order.EQ.5) then
  
    Sfilter0%b%coeff(1)=3.236d0 
    Sfilter0%b%coeff(2)=5.236d0 
    Sfilter0%b%coeff(3)=5.236d0 
    Sfilter0%b%coeff(4)=3.236d0 
 
  else if (order.EQ.6) then

    Sfilter0%b%coeff(1)=3.864d0 
    Sfilter0%b%coeff(2)=7.464d0 
    Sfilter0%b%coeff(3)=9.142d0 
    Sfilter0%b%coeff(4)=7.464d0 
    Sfilter0%b%coeff(5)=3.864d0 
 
  else if (order.EQ.7) then

    Sfilter0%b%coeff(1)=4.494d0 
    Sfilter0%b%coeff(2)=10.098d0 
    Sfilter0%b%coeff(3)=14.592d0 
    Sfilter0%b%coeff(4)=14.592d0 
    Sfilter0%b%coeff(5)=10.098d0 
    Sfilter0%b%coeff(6)=4.494d0 
 
  else if (order.EQ.8) then

    Sfilter0%b%coeff(1)=5.126d0 
    Sfilter0%b%coeff(2)=13.137d0 
    Sfilter0%b%coeff(3)=21.846d0 
    Sfilter0%b%coeff(4)=25.688d0 
    Sfilter0%b%coeff(5)=21.846d0 
    Sfilter0%b%coeff(6)=13.137d0 
    Sfilter0%b%coeff(7)=5.126d0 
 
  else if (order.EQ.9) then

    Sfilter0%b%coeff(1)=5.759d0 
    Sfilter0%b%coeff(2)=16.582d0 
    Sfilter0%b%coeff(3)=31.163d0 
    Sfilter0%b%coeff(4)=41.986d0 
    Sfilter0%b%coeff(5)=41.986d0 
    Sfilter0%b%coeff(6)=31.163d0 
    Sfilter0%b%coeff(7)=16.582d0 
    Sfilter0%b%coeff(8)=5.759d0 
 
  else if (order.EQ.10) then

    Sfilter0%b%coeff(1)=6.392d0 
    Sfilter0%b%coeff(2)=20.432d0 
    Sfilter0%b%coeff(3)=42.802d0 
    Sfilter0%b%coeff(4)=64.882d0 
    Sfilter0%b%coeff(5)=74.233d0 
    Sfilter0%b%coeff(6)=64.882d0 
    Sfilter0%b%coeff(7)=42.802d0 
    Sfilter0%b%coeff(8)=20.432d0 
    Sfilter0%b%coeff(9)=6.392d0 

  end if

  write(record_user_inputs_unit,'(A)')ch
  write(post_process_info_unit,'(A13,A,A2)')'Filter type: ',ch,'pf'
  
  write(record_user_inputs_unit,'(ES16.6)')frequency
  write(post_process_info_unit,'(A18,ES16.6)')'Cutoff frequency: ',frequency
  
  write(record_user_inputs_unit,'(I6)')order
  write(post_process_info_unit,'(A14,I6)')'Filter order: ',order

! Write the filter function to file

  write(*,*)'Enter the filename for the filter function'
  read(*,'(A)')filter_file_name

  open(UNIT=local_file_unit, 					    &
       FILE=trim(filter_file_name),	    &
       ERR=9000)

  call write_Sfilter(Sfilter0,local_file_unit) ! filter function

  close(UNIT=local_file_unit)
  
  write(record_user_inputs_unit,'(A)')trim(filter_file_name)  
  write(post_process_info_unit,*)'	Filter filename:',trim(filter_file_name)

! deallocate memory  
  CALL deallocate_Sfilter( Sfilter0 )
  
  RETURN
  
9000 CALL write_line('Error writing filter function',0,.TRUE.)
     CALL write_line('Problem opening filter file:'//trim(filter_file_name),0,.TRUE.)
     STOP
  
END SUBROUTINE create_filter_function
