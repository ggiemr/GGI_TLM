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
! SUBROUTINE plot_time_domain_data
!
! NAME
!    plot_time_domain_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 17/08/2012 CJS
!
!
SUBROUTINE plot_time_domain_data

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: data_filename
  character(len=256)	:: gnuplot_filename
  character(len=256)	:: curve_title
  character(len=256)	:: y_label
  character		:: plot_option
  
  integer		:: timestep
  
  integer		:: function_number
  
  logical		:: file_exists

  character(len=256)	:: command

! START

  n_functions_of_time=0
  n_functions_of_frequency=0

5 write(*,*)
  write(*,*)'Time domain output files:'
  
  command='ls -ltr *.tout'
  CALL system(command)

  write(*,*)'Enter the time domain filename'
  read(*,*)data_filename
  inquire(file=trim(data_filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(data_filename)

  write(*,*)'Enter the label for the y axis'
  read(*,'(A)')y_label
  write(record_user_inputs_unit,'(A)')trim(y_label)
  
  write(*,*)'Enter the title for the time domain curve'
  read(*,'(A)')curve_title
  write(record_user_inputs_unit,'(A)')trim(curve_title)
  
  write(*,*)'Enter the output required s(creen) g(if) j(peg)'
  read(*,'(A)')plot_option
  write(record_user_inputs_unit,'(A)')trim(plot_option)
  
  CALL convert_to_lower_case(plot_option,1)
    
! Write gnuplot file to plot this data    

  gnuplot_filename='Time_domain_plot.plt'
  OPEN(unit=local_file_unit,file=gnuplot_filename)
  
  if (plot_option.eq.'g') then
    write(local_file_unit,'(A)')'set term gif'
    write(local_file_unit,'(A)')'set output "time_domain_data_plot.gif"'   
  else if (plot_option.eq.'j') then
    write(local_file_unit,'(A)')'set term jpeg'
    write(local_file_unit,'(A)')'set output "time_domain_data_plot.jpeg"'   
  end if
  
  write(local_file_unit,'(A)')'set xlabel "Time (s)"'
  write(local_file_unit,'(A,A,A)')'set ylabel "',trim(y_label),'"'
  write(local_file_unit,'(A,A,A,A,A)')'plot "',trim(data_filename),'" u 1:3 title"',trim(curve_title),'" w l'
  
  if (plot_option.eq.'s') then
    write(local_file_unit,'(A)')'pause 100'
  end if
  
  CLOSE(unit=local_file_unit)
    
  command='gnuplot Time_domain_plot.plt & '
  CALL system(command)

  RETURN
  
END SUBROUTINE plot_time_domain_data
