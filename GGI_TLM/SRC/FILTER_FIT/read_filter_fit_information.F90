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
! SUBROUTINE read_filter_fit_information
!
! NAME
!     read filter fit information
!
! DESCRIPTION
!     read the options and data associated with the filter fitting process
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/12/2012 CJS
!
!
SUBROUTINE read_filter_fit_information

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff

IMPLICIT NONE
 
! local variables

  character	:: ch_in

! START

  write(*,*)'Enter the name of the filter fitting input file'
  read(*,'(A256)')FF_name
 
  CALL write_line('Problem name:'//trim(FF_name),0,ff_output_to_screen)

  write(*,*)'Enter model type:'
  write(*,*)'1 : dielectric material model fit '
  write(*,*)'2 : magnetic material model fit '
  write(*,*)'3 : thin layer model fit'
  write(*,*)'4 : impedance model fit '
  
  read(*,*)fit_type
  
  if ( (fit_type.lt.1).or.(fit_type.gt.4) ) then
    write(*,*)'Error specifying fit_type'
    write(*,*)'Fit_type=',fit_type
    write(*,*)   &
    'Fittype must be 1 to 4'
    STOP
  end if
         
! read frequency scaling factor      
  write(*,*)'Enter frequency scaling factor to Hz'
  read(*,*)fscale
  write(*,*)'Frequency scaling factor=',fscale
       
! read data
  call read_data()
       
  write(*,*)'Enter model order to fit'
  read(*,*)order
  
  stabilise_input_data_flag=.FALSE.
  write(*,*)'Do you want to stabilise the input data (y/n)'
  read(*,'(A1)')ch_in
  if ( (ch_in.eq.'y').OR.(ch_in.eq.'Y') ) then
    stabilise_input_data_flag=.TRUE.
  end if
  
  stabilise_filter_flag=.FALSE.
  write(*,*)'Do you want to stabilise the filter before optimisation (y/n)'
  read(*,'(A1)')ch_in
  if ( (ch_in.eq.'y').OR.(ch_in.eq.'Y') ) then
    stabilise_filter_flag=.TRUE.
  end if
  
  optimise_filter_flag=.FALSE.
  write(*,*)'Do you want to optimise the filter function (y/n)'
  read(*,'(A1)')ch_in
  if ( (ch_in.eq.'y').OR.(ch_in.eq.'Y') ) then
    optimise_filter_flag=.TRUE.
  end if

  RETURN

END SUBROUTINE read_filter_fit_information

