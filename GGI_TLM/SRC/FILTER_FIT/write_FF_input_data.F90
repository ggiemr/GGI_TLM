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
! SUBROUTINE write_FF_input_data
!
! NAME
!     write_FF_input_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/12/2012 CJS
!
!
SUBROUTINE write_FF_input_data

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff

IMPLICIT NONE
 
! local variables

  integer	:: function_loop
  integer	:: freq_loop
  
  character(len=256) :: filename(1:4)

! START

  CALL write_line('CALLED: write_FF_input_data',0,ff_output_to_screen)

  if (fit_type.eq.dielectric_material) then 
  
    filename(1)=trim(FF_name)//dielectric_input_data_extension
 
  else if (fit_type.eq.magnetic_material) then 
  
    filename(1)=trim(FF_name)//magnetic_input_data_extension
 
  else if (fit_type.eq.thin_layer) then 
  
    filename(1)=trim(FF_name)//thin_layer_z11_input_data_extension
    filename(2)=trim(FF_name)//thin_layer_z12_input_data_extension
    filename(3)=trim(FF_name)//thin_layer_z21_input_data_extension
    filename(4)=trim(FF_name)//thin_layer_z22_input_data_extension
 
  else if (fit_type.eq.impedance) then 
   
    filename(1)=trim(FF_name)//impedance_input_data_extension

  end if
  
  
  do function_loop=1,n_functions
  
    OPEN(unit=input_data_file_unit,file=filename(function_loop))
  
    do freq_loop=1,n_values

      write(input_data_file_unit,8000)frequency(freq_loop)  ,	&
                                      dble(value(function_loop,freq_loop)),	&
                                      dimag(value(function_loop,freq_loop))
8000  format(3E16.6)
     
    end do ! next freq_loop
    
  end do ! next function

  CALL write_line('FINISHED: write_FF_input_data',0,ff_output_to_screen)

END SUBROUTINE write_FF_input_data
