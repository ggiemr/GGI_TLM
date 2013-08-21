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
! SUBROUTINE read_SAR_volume_list
!
! NAME
!     read_SAR_volume_list
!
! DESCRIPTION
!     read SAR volume list packet
!
! Example packet:
!                        
!SAR_volume_list
!2 ! number of SAR volumes
!1  ! SAR VOLUME NUMBER
!1  ! material number
!900E6  ! frequency
!968.6  ! material density kg/m^3
!2  ! SAR VOLUME NUMBER
!2  ! material number
!900E6  ! frequency
!968.6  ! material density kg/m^3
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 9/01/2013 CJS
!
!
SUBROUTINE read_SAR_volume_list

USE TLM_general
USE file_information
USE geometry
USE TLM_output
USE cell_parameters

IMPLICIT NONE

! local variables

integer	:: output_number
integer	:: read_number

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_SAR_volume_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_SAR_volumes
  
  CALL write_line_integer('number of output surfaces',n_SAR_volumes,0,output_to_screen_flag)
  
  if ( allocated( SAR_volume_list ) ) GOTO 9000
  
  ALLOCATE ( SAR_volume_list(1:n_SAR_volumes) )

  do output_number=1,n_SAR_volumes
  
    CALL write_line_integer('Reading output number',output_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.output_number) goto 9010
    
    read(input_file_unit,*,err=9005)SAR_volume_list(output_number)%material_number
    
    read(input_file_unit,*,err=9005)SAR_volume_list(output_number)%frequency
    
    read(input_file_unit,*,err=9005)SAR_volume_list(output_number)%density
    
  end do ! next SAR volume

  CALL write_line('FINISHED: read_SAR_volume_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating SAR_volume_list:',0,.TRUE.)
     CALL write_line('SAR_volume_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading SAR volume list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading SAR volume list packet',0,.TRUE.)
     CALL write_line('SAR volumes should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_SAR_volume_list
