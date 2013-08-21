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
! SUBROUTINE write_time_domain_header_data
! SUBROUTINE read_time_domain_header_data
! SUBROUTINE write_frequency_domain_header_data
! SUBROUTINE read_frequency_domain_header_data
!
! NAME
!     write_time_domain_header_data
!
! DESCRIPTION
!     write header information relevant to time domain output files
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/08/2012 CJS
!
!
SUBROUTINE write_time_domain_header_data(file_unit,n_data_points,n_timesteps)

USE file_header

IMPLICIT NONE

  integer 	:: file_unit
  integer 	:: n_data_points
  integer 	:: n_timesteps

! local variables

! START
  
!  CALL write_line('CALLED: write_time_domain_header_data',0,output_to_screen_flag)

   write(file_unit,'(A)')time_domain_output_file_header
   write(file_unit,'(A25,I10)')  '# Number of data points: ',n_data_points
   write(file_unit,'(A25,I10)')  '# Number of Timesteps  : ',n_timesteps
  
!  CALL write_line('FINISHED: write_time_domain_header_data',0,output_to_screen_flag)

  RETURN

END SUBROUTINE write_time_domain_header_data
!
! NAME
!     read_time_domain_header_data
!
! DESCRIPTION
!     read header information relevant to time domain output files
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/08/2012 CJS
!
!
SUBROUTINE read_time_domain_header_data(file_unit,n_data_points,n_timesteps)

USE file_header

IMPLICIT NONE

  integer 	:: file_unit
  integer 	:: n_data_points
  integer 	:: n_timesteps

! local variables

  character(len=34)	:: title
  character(len=25)	:: description

! START
  
!  CALL write_line('CALLED: read_time_domain_header_data',0,output_to_screen_flag)

   read(file_unit,'(A34)')title
   
! check title string   
   if (title.ne.time_domain_output_file_header) goto 9000
   
   read(file_unit,'(A25,I10)')  description,n_data_points
   read(file_unit,'(A25,I10)')  description,n_timesteps
  
!  CALL write_line('FINISHED: read_time_domain_header_data',0,output_to_screen_flag)

  RETURN

9000  CALL write_line('Error in read_time_domain_header_data',0,.TRUE.)
      CALL write_line('Not a time domain file',0,.TRUE.)
      CALL write_line('Expecting a header title:',0,.TRUE.)
      CALL write_line(time_domain_output_file_header,0,.TRUE.)
      STOP
 
END SUBROUTINE read_time_domain_header_data
!
! NAME
!     write_frequency_domain_header_data
!
! DESCRIPTION
!     write header information relevant to frequency domain output files
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/08/2012 CJS
!
!
SUBROUTINE write_frequency_domain_header_data(file_unit,n_data_points,n_frequencies)

USE file_header

IMPLICIT NONE

  integer 	:: file_unit
  integer 	:: n_data_points
  integer 	:: n_frequencies

! local variables

! START
  
!  CALL write_line('CALLED: write_frequency_domain_header_data',0,output_to_screen_flag)

   write(file_unit,'(A)')frequency_domain_output_file_header
   write(file_unit,'(A25,I10)')  '# Number of data points: ',n_data_points
   write(file_unit,'(A25,I10)')  '# Number of frequencies: ',n_frequencies
  
!  CALL write_line('FINISHED: write_frequency_domain_header_data',0,output_to_screen_flag)

  RETURN

END SUBROUTINE write_frequency_domain_header_data
!
! NAME
!     read_frequency_domain_header_data
!
! DESCRIPTION
!     read header information relevant to frequency domain output files
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/08/2012 CJS
!
!
SUBROUTINE read_frequency_domain_header_data(file_unit,n_data_points,n_frequencies)

USE file_header

IMPLICIT NONE

  integer 	:: file_unit
  integer 	:: n_data_points
  integer 	:: n_frequencies

! local variables

  character(len=34)	:: title
  character(len=25)	:: description

! START
  
!  CALL write_line('CALLED: read_frequency_domain_header_data',0,output_to_screen_flag)

   read(file_unit,'(A34)')title
   
! check title string   
   if (title.ne.frequency_domain_output_file_header) goto 9000
   
   read(file_unit,'(A25,I10)')  description,n_data_points
   read(file_unit,'(A25,I10)')  description,n_frequencies
  
!  CALL write_line('FINISHED: read_frequency_domain_header_data',0,output_to_screen_flag)

  RETURN

9000  CALL write_line('Error in read_frequency_domain_header_data',0,.TRUE.)
      CALL write_line('Not a frequency domain file',0,.TRUE.)
      CALL write_line('Expecting a header title:',0,.TRUE.)
      CALL write_line(frequency_domain_output_file_header,0,.TRUE.)
      STOP
 
END SUBROUTINE read_frequency_domain_header_data
