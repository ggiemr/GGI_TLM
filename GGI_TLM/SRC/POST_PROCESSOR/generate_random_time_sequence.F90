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
! SUBROUTINE Generate_random_sequence
!
! NAME
!    Generate_random_sequence
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/11/2013 CJS
!
!
SUBROUTINE Generate_random_sequence()

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  integer		:: n_sequences
  integer		:: n_timesteps
  
  real*8		:: dt
  
  integer		:: function_number
  integer		:: loop
  integer		:: timestep
  real*8 		:: r1,r2,value,value2
  real*8		:: mean
  real*8		:: sigma
  real*8		:: range,xmin,xmax
  
  integer		:: distribution_type
  integer,parameter	:: uniform=1
  integer,parameter	:: normal=2
  
  logical		:: normal_value_calculated
  
  character 	:: ch
  logical	:: correlation_flag
  
  integer	:: row,col
  
  real*8,allocatable		:: C(:,:)
  real*8,allocatable		:: L(:,:)
  real*8,allocatable		:: U(:,:)

  real*8,allocatable		:: F1(:)
  real*8,allocatable		:: F2(:)
  
! START

  write(*,*)
  write(*,*)'Generate_random_sequence'
  
  write(*,*)'Enter the number of random sequences to generate'
  read(*,*)n_sequences
  write(record_user_inputs_unit,*)n_sequences,' Number of random sequences'
  
  n_functions_of_time=n_sequences
  n_functions_of_frequency=0
  
  CALL Allocate_post_data()
  
  write(*,*)'Enter the number of timesteps to generate'
  read(*,*)n_timesteps
  write(record_user_inputs_unit,*)n_timesteps,' Number of timesteps'

  write(*,*)'Enter the timestep'
  read(*,*)dt
  write(record_user_inputs_unit,*)dt,' timestep'
  
  write(*,*)'Distribution from which the samples are to be drawn: '
  write(*,*)'Uniform or Normal (u or n)'
  read(*,'(A)')ch
  if ( (ch.eq.'u').OR.(ch.eq.'U') ) then
    distribution_type=uniform
  else if ( (ch.eq.'n').OR.(ch.eq.'N') ) then
    distribution_type=normal
  else
    write(*,*)"Response should be 'u' or 'n'"
    STOP
  end if
  write(record_user_inputs_unit,'(A)')ch
  
  if (distribution_type.EQ.uniform)then
  
    write(*,*)'Enter the minimum value for the uniformly distributed values'
    read(*,*)xmin
    write(record_user_inputs_unit,*)xmin,' minimum value for the uniformly distributed values'
  
    write(*,*)'Enter the maximum value for the uniformly distributed values'
    read(*,*)xmax
    write(record_user_inputs_unit,*)xmax,' maximum value for the uniformly distributed values'
    
    mean=(xmax+xmin)/2d0
    range=xmax-xmin
    
  else if (distribution_type.EQ.normal)then
  
    write(*,*)'Enter the mean value'
    read(*,*)mean
    write(record_user_inputs_unit,*)mean,' mean'
  
    write(*,*)'Enter the standard deviation, sigma'
    read(*,*)sigma
    write(record_user_inputs_unit,*)sigma,' sigma'
    
  end if
  
  normal_value_calculated=.FALSE.
  
  do function_number=1,n_sequences
  
    function_of_time(function_number)%n_timesteps=n_timesteps
    ALLOCATE ( function_of_time(function_number)%time(1:n_timesteps) )
    ALLOCATE ( function_of_time(function_number)%value(1:n_timesteps) )
    
    do timestep=1,n_timesteps
      function_of_time(function_number)%time(timestep)=dt*(timestep-1)
      
      if (distribution_type.EQ.uniform) then
      
        CALL random_number(r1)
	value=(r1-0.5d0)*range+mean
	
      else ! normal distribution
      
        if (.NOT.normal_value_calculated) then
! calculates two independent normally distributed samples at a time
          CALL random_number(r1)
	  CALL random_number(r2)
	  value =sqrt(-2d0*log(r1))*cos(6.28318530718d0*r2)
	  value2=sqrt(-2d0*log(r1))*sin(6.28318530718d0*r2)
	  normal_value_calculated=.TRUE.
	else
	  value=value2
	  normal_value_calculated=.FALSE.
	end if
	
	value=value*sigma+mean
	
      end if
      	
      function_of_time(function_number)%value(timestep)=value
      
    end do
    
  end do
  
110 CONTINUE
  write(*,*)'Correlate source sequences? (y or n)'
  read(*,'(A)')ch
  if ( (ch.eq.'y').OR.(ch.eq.'Y') ) then
    correlation_flag=.TRUE.
  else if ( (ch.eq.'n').OR.(ch.eq.'N') ) then
    correlation_flag=.FALSE.
  else
    write(*,*)"Response should be 'y' or 'n'"
    GOTO 110
  end if
  write(record_user_inputs_unit,'(A)')ch
  
  if (correlation_flag) then
  
    ALLOCATE( C(1:n_sequences,1:n_sequences) )
    ALLOCATE( L(1:n_sequences,1:n_sequences) )
    ALLOCATE( U(1:n_sequences,1:n_sequences) )
    ALLOCATE( F1(1:n_sequences) )
    ALLOCATE( F2(1:n_sequences) )
  
    write(*,*)'Please enter the correlation matrix'
  
    do row=1,n_sequences
    
      do col=row+1,n_sequences
  
        write(*,*)'Enter matrix element ',row,col
	read(*,*)C(row,col)
        write(record_user_inputs_unit,*)C(row,col)
	C(col,row)=C(row,col)
	
      end do
      
      C(row,row)=1d0
      
    end do

    CALL Cholesky(C,n_sequences,L,U)
    
    do timestep=1,n_timesteps

      do row=1,n_sequences
        F1(row)=function_of_time(row)%value(timestep)
      end do

      CALL dmatvmul(L,n_sequences,n_sequences,F1,n_sequences,F2,n_sequences)
      
      do row=1,n_sequences
        function_of_time(row)%value(timestep)=F2(row)
      end do
      
    end do

    DEALLOCATE( C )
    DEALLOCATE( L )
    DEALLOCATE( U )
    DEALLOCATE( F1 )
    DEALLOCATE( F2 )
    
  end if ! correlation_flag
	    
  CALL write_time_domain_data_array(n_sequences)
  
! Deallocate memory  
 
  CALL Deallocate_post_data()
  
  RETURN
  
  
END SUBROUTINE Generate_random_sequence
