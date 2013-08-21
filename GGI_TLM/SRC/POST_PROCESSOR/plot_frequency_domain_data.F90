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
! SUBROUTINE plot_frequency_domain_data
!
! NAME
!    plot_frequency_domain_data
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
SUBROUTINE plot_frequency_domain_data


USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: data_filename
  character(len=256)	:: data_filename_with_quotes
  character(len=256)	:: filename
  character(len=256)	:: curve_title
  character(len=256)	:: y_label
  character		:: plot_option
  character		:: ch
  
  integer		:: frequency
  
  logical		:: plot_real
  logical		:: plot_imag
  logical		:: plot_mag
  logical		:: plot_phase
  logical		:: plot_dB
  
  logical		:: plot_linear
  
  character		:: first_character
  
  integer		:: function_number
  
  logical		:: file_exists

  character(len=256)	:: command
  
! START

5 write(*,*)
  write(*,*)'Frequency domain output files:'
  
  command='ls -ltr *.fout'
  CALL system(command)

  write(*,*)'Enter the frequency domain filename'
  read(*,*)data_filename
  inquire(file=trim(data_filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(data_filename)
  
  plot_linear=.FALSE.
  plot_dB=.FALSE.

  write(*,*)'plot real part? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  plot_real=.FALSE.
  if (ch.eq.'y') then
    plot_real=.TRUE.
    plot_linear=.TRUE.
  end if
  
  write(*,*)'plot imaginary part? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  plot_imag=.FALSE.
  if (ch.eq.'y') then
    plot_imag=.TRUE.
    plot_linear=.TRUE.
  end if

  write(*,*)'plot magnitude? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  plot_mag=.FALSE.
  if (ch.eq.'y') then
    plot_mag=.TRUE.
    plot_linear=.TRUE.
  end if

  write(*,*)'plot phase? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  plot_phase=.FALSE.
  if (ch.eq.'y') then
    plot_phase=.TRUE.
    plot_linear=.TRUE.
  end if

  write(*,*)'plot magnitude in dB? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  if (ch.eq.'y') then
    plot_dB=.TRUE.
  end if

  write(*,*)'Enter the label for the y axis'
  read(*,'(A)')y_label
  write(record_user_inputs_unit,'(A)')trim(y_label)
  
  write(*,*)'Enter the title for the frequency domain domain curves'
  read(*,'(A)')curve_title
  write(record_user_inputs_unit,'(A)')trim(curve_title)
  
  write(*,*)'Enter the output required s(creen) g(if) j(peg)'
  read(*,'(A)')plot_option
  write(record_user_inputs_unit,'(A)')trim(plot_option)
  
  CALL convert_to_lower_case(plot_option,1)
  
  data_filename_with_quotes='"'//trim(data_filename)//'"'

  if (plot_linear) then

    filename='Frequency_domain_plot_lin.plt'
    OPEN(unit=local_file_unit,file=filename)
  
    if (plot_option.eq.'g') then
      write(local_file_unit,'(A)')'set term gif'
      write(local_file_unit,'(A)')'set output "frequency_domain_data_plot.gif"'   
    else if (plot_option.eq.'j') then
      write(local_file_unit,'(A)')'set term jpeg'
      write(local_file_unit,'(A)')'set output "frequency_domain_data_plot.jpeg"'   
    end if
  
    if (plot_phase) then
  
      write(local_file_unit,'(A)')'set ytics nomirror'
      write(local_file_unit,'(A)')'set y2range[-3.142:3.142]'
      write(local_file_unit,'(A)')'set y2tics'
      write(local_file_unit,'(A)')'show y2tics'
      write(local_file_unit,'(A)')'set y2label "Phase (radians)"'

    end if
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A,A,A)')'set ylabel "',trim(y_label),'"'
    write(local_file_unit,'(A)')'plot \'
  
    first_character=' '
  
    if (plot_real) then
      write(local_file_unit,'(A,A,A,A,A)')first_character,trim(data_filename_with_quotes),' u 1:3 title"Re{', &
                                                        trim(curve_title),'}" w l \'
      first_character=','
    end if
    if (plot_imag) then
      write(local_file_unit,'(A,A,A,A,A)')first_character,trim(data_filename_with_quotes),' u 1:4 title"Im{', &
                                                        trim(curve_title),'}" w l \'
      first_character=','
    end if
    if (plot_mag) then
      write(local_file_unit,'(A,A,A,A,A)')first_character,trim(data_filename_with_quotes),' u 1:5 title"|', &
                                                        trim(curve_title),'|" w l \'
      first_character=','
    end if
    if (plot_phase) then
      write(local_file_unit,'(A,A,A,A,A)')first_character,trim(data_filename_with_quotes),' u 1:6 axes x1y2 title"phase (', &
                                                        trim(curve_title),')" w l \'
    end if
  
    write(local_file_unit,'(A)')' '
   
    if (plot_option.eq.'s') then
      write(local_file_unit,'(A)')'pause 100'
    end if
  
    CLOSE(unit=local_file_unit)
    
    command='gnuplot Frequency_domain_plot_lin.plt & '
    CALL system(command)

  end if
  
  if (plot_dB) then

    filename='Frequency_domain_plot_dB.plt'
    OPEN(unit=local_file_unit,file=filename)
  
    if (plot_option.eq.'g') then
      write(local_file_unit,'(A)')'set term gif'
      write(local_file_unit,'(A)')'set output "frequency_domain_data_plot_dB.gif"'   
    else if (plot_option.eq.'j') then
      write(local_file_unit,'(A)')'set term jpeg'
      write(local_file_unit,'(A)')'set output "frequency_domain_data_plot_dB.jpeg"'   
    end if
    
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A,A,A)')'set ylabel "',trim(y_label),' (dB)"'
    write(local_file_unit,'(A)')'plot \'
  
    first_character=' '
  
    write(local_file_unit,'(A,A,A,A,A)')first_character,trim(data_filename_with_quotes),' u 1:7 title"', &
                                                        trim(curve_title),'" w l '
  
    write(local_file_unit,'(A)')' '
   
    if (plot_option.eq.'s') then
      write(local_file_unit,'(A)')'pause 100'
    end if
  
    CLOSE(unit=local_file_unit)
    
    command='gnuplot Frequency_domain_plot_dB.plt & '
    CALL system(command)

  end if
   

  RETURN
  

  
END SUBROUTINE plot_frequency_domain_data
