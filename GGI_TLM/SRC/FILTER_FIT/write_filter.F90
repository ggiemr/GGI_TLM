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
! SUBROUTINE write_filter
!
! NAME
!     write_filter
!
! DESCRIPTION
!     write filters to material files in the appropriate format
!     
! COMMENTS
!     include a header string:
!# NAME, Model order= order, stable/unstable
!
! HISTORY
!
!     started 11/12/2012 CJS
!
!
SUBROUTINE write_filter

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff

USE filter_types
USE filter_functions

USE constants

IMPLICIT NONE
 
! local variables

  character(len=256) :: temp_filename
  character(len=256) :: filename

  character(len=256) :: header_string
  character(len=256) :: temp_string
  character(len=256) :: stability_string
  character(len=256) :: model_order_string
  
  type(Sfilter)		:: unit_filter

! START

  CALL write_line('CALLED: write_filter_frequency_response',0,ff_output_to_screen)

! create a unit filter
  unit_filter=allocate_Sfilter(0,0)
  unit_filter%a%coeff(0)=1d0
  unit_filter%b%coeff(0)=1d0

! get the material filter filename  
  CALL add_integer_to_filename(FF_name,order,temp_filename)

  if (fit_type.eq.dielectric_material) then 
    filename=trim(temp_filename)//dielectric_filter_extension
  else if (fit_type.eq.magnetic_material) then  
    filename=trim(temp_filename)//magnetic_filter_extension
  else if (fit_type.eq.thin_layer) then  
    filename=trim(temp_filename)//thin_layer_filter_extension
  else if (fit_type.eq.impedance) then  
    filename=trim(temp_filename)//impedance_filter_extension
  end if

! create the header string  
  if (stable_filter) then
    stability_string=' : STABLE'
  else
    stability_string=' : UNSTABLE'
  end if 
  
  temp_string='# '//trim(FF_name)//' Model order '
   
  CALL add_integer_to_filename(temp_string,order,model_order_string)  
  
  header_string=trim(model_order_string)//'. created by GGI_filter_fit '//trim(stability_string)

! open the material filter file  
  OPEN(unit=filter_file_unit,file=filename)
  
  write(filter_file_unit,'(A)')trim(header_string)

  write(filter_file_unit,*)frequency(1),frequency(n_values),' # frequency range of validity'
  
  if (fit_type.eq.dielectric_material) then 
    
    write(filter_file_unit,'(A)')'# permittivity filter'
    CALL write_Sfilter(filter_S(1),filter_file_unit)
    write(filter_file_unit,*)filter_sigma(1)*wnorm*eps0,' # electric conductivity'
    
    write(filter_file_unit,'(A)')'# permeability filter'
    CALL write_Sfilter(unit_filter,filter_file_unit)
    write(filter_file_unit,*)0d0,' # magnetic conductivity'
    
  else if (fit_type.eq.magnetic_material) then  
    
    write(filter_file_unit,'(A)')'# permittivity filter'
    CALL write_Sfilter(unit_filter,filter_file_unit)
    write(filter_file_unit,*)0d0,' # electric conductivity'
    
    write(filter_file_unit,'(A)')'# permeability filter'
    CALL write_Sfilter(filter_S(1),filter_file_unit)
    write(filter_file_unit,*)filter_sigma(1)*wnorm*mu0,' # magnetic conductivity'
    
  else if (fit_type.eq.thin_layer) then  
    
    write(filter_file_unit,'(A)')'# z11 filter'
    CALL write_Sfilter(filter_S(1),filter_file_unit)   
    write(filter_file_unit,'(A)')'# z12 filter'
    CALL write_Sfilter(filter_S(2),filter_file_unit)    
    write(filter_file_unit,'(A)')'# z21 filter'
    CALL write_Sfilter(filter_S(3),filter_file_unit)    
    write(filter_file_unit,'(A)')'# z22 filter'
    CALL write_Sfilter(filter_S(4),filter_file_unit)
       
  else if (fit_type.eq.impedance) then  
    
    write(filter_file_unit,'(A)')'# Impedance filter'
    CALL write_Sfilter(filter_S(1),filter_file_unit)
    
  end if
  
  CLOSE(unit=filter_file_unit)
  
  CALL deallocate_Sfilter( unit_filter )
  
  CALL write_line('FINISHED: write_filter_frequency_response',0,ff_output_to_screen)

  RETURN

END SUBROUTINE write_filter
