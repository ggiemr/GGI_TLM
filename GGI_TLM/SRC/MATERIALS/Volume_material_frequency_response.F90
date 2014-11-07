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
!
! Name Volume_material_frequency_response
!     
!
! Description:
!      Calculate the frequency response of volume materials
!      
!
! Comments:
!      
!      
!
! History
!
!     started 07/11/2012 CJS
!     17/10/2014 CJS - add write_all_material_info_to_file flag - eliminates the need for user input in this case
!

SUBROUTINE Volume_material_frequency_response(write_all_material_info_to_file)

USE TLM_general
USE TLM_volume_materials
USE filter_types
USE filter_operators
USE filter_functions
USE file_general
USE file_information
USE constants

IMPLICIT NONE

! variables passed to subroutine
  
  logical :: write_all_material_info_to_file

! local_variables
  
  integer :: material_number,first_material,last_material
  character(LEN=filename_length)    :: eps_filename
  character(LEN=filename_length)    :: mu_filename
  real*8  :: fmin,fmax,fstep
  real*8  :: sigma
  character(LEN=filename_length)    :: gnuplot_filename
  character(LEN=filename_length)    :: base_filename
  character(LEN=filename_length)    :: jpeg_filename
  character(LEN=line_length)	    :: command

! function_types

! START

  write(*,*)'CALLED: Volume_material_frequency_response'
  write(*,*)
  write(*,*)'Number of volume materials=',n_volume_materials
  
  if (n_volume_materials.EQ.0) then
    write(*,*)'No volume materials to view'
    RETURN
  end if

! read the material number to output

! read the volume material number to view  
  if (write_all_material_info_to_file) then
! set the material number to zero to indicate that all material data should be written
    material_number=0
  else
    write(*,*)
    write(*,*)'Enter the volume material number to view or 0 to view all of them'
    read(*,*)material_number
  end if
  
  if ( (material_number.lt.0).OR.(material_number.gt.n_volume_materials) ) then
    write(*,*)'volume_material_number is outside the available range'
    write(*,*)'Number of volume materials=',n_volume_materials
    RETURN
  end if
  
  if (material_number.eq.0) then
    first_material=1
    last_material=n_volume_materials
  else
    first_material=material_number
    last_material=material_number
  end if

  do material_number=first_material,last_material

    if (volume_material_list(material_number)%type.EQ.volume_material_type_PEC) then
  
      write(*,*)'Material number',material_number,' is PEC'
  
    else if (volume_material_list(material_number)%type.EQ.volume_material_type_PMC) then
  
      write(*,*)'Material number',material_number,' is PMC'
  
    else if (volume_material_list(material_number)%type.EQ.volume_material_type_DISPERSIVE) then
    
      fmin=volume_material_list(material_number)%fmin
      fmax=volume_material_list(material_number)%fmax
      fstep=(fmax-fmin)/200d0

! Write Dielectric filter response
      eps_filename=trim(volume_material_list(material_number)%name)//'.eps.fout'

! get the jpg output filenames - if we are running automatically then tag with the material number
      if (write_all_material_info_to_file) then
        base_filename=trim(volume_material_list(material_number)%name)//'_eps.jpg'
        CALL add_integer_to_filename(base_filename,material_number,jpeg_filename)
      else
        jpeg_filename=trim(volume_material_list(material_number)%name)//'_eps.jpg'
      end if
      
      sigma=volume_material_list(material_number)%sigma_e
    
      CALL output_material_frequency_response(volume_material_list(material_number)%eps_S,	&
                                              sigma,fmin,fmax,fstep,local_file_unit,eps_filename)
      
! Write gnuplot file

      gnuplot_filename=trim(volume_material_list(material_number)%name)//'_eps.plt'
      OPEN(unit=local_file_unit,file=gnuplot_filename)
  
      write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
      write(local_file_unit,'(A)')'set ylabel "Relative Permittivity"'
      write(local_file_unit,'(A,A,A)')'set title "',trim(volume_material_list(material_number)%name),'"'

      write(local_file_unit,'(A)')'set term jpeg'
      write(local_file_unit,'(A,A,A)')'set output "',trim(jpeg_filename),'"'
      write(local_file_unit,'(A,A,A)')'plot "',trim(eps_filename),'" u 1:2 title "Re{epsr}" w lp,\'
      write(local_file_unit,'(A,A,A)')'     "',trim(eps_filename),'" u 1:3 title "Im{epsr}" w lp'

      if (.NOT.write_all_material_info_to_file) then
! only plot to the screen if we are operating in the interactive mode
        write(local_file_unit,'(A)')'set term wxt 0'
        write(local_file_unit,'(A,A,A)')'plot "',trim(eps_filename),'" u 1:2 title "Re{epsr}" w lp,\'
        write(local_file_unit,'(A,A,A)')'     "',trim(eps_filename),'" u 1:3 title "Im{epsr}" w lp'
        write(local_file_unit,'(A)')'pause 100'
      end if
      
      CLOSE(unit=local_file_unit)

! plot permittivity response

      command='gnuplot '//trim(gnuplot_filename)//' & '
      CALL system(command)
      
! Write Magnetic filter response
  
      mu_filename=trim(volume_material_list(material_number)%name)//'.mu.fout'

! get the jpg output filenames - if we are running automatically then tag with the material number
      if (write_all_material_info_to_file) then
        base_filename=trim(volume_material_list(material_number)%name)//'_mu.jpg'
        CALL add_integer_to_filename(base_filename,material_number,jpeg_filename)
      else
        jpeg_filename=trim(volume_material_list(material_number)%name)//'_mu.jpg'
      end if
    
      sigma=volume_material_list(material_number)%sigma_e
    
      CALL output_material_frequency_response(volume_material_list(material_number)%mu_S,	&
                                              sigma,fmin,fmax,fstep,local_file_unit,mu_filename)

! Write gnuplot file

      gnuplot_filename=trim(volume_material_list(material_number)%name)//'_mu.plt'
      OPEN(unit=local_file_unit,file=gnuplot_filename)
  
      write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
      write(local_file_unit,'(A)')'set ylabel "Relative Permeability"'
      write(local_file_unit,'(A,A,A)')'set title "',trim(volume_material_list(material_number)%name),'"'

      write(local_file_unit,'(A)')'set term jpeg'
      write(local_file_unit,'(A,A,A)')'set output "',trim(jpeg_filename),'"'
      write(local_file_unit,'(A,A,A)')'plot "',trim(mu_filename),'" u 1:2 title "Re{mur}" w lp,\'
      write(local_file_unit,'(A,A,A)')'     "',trim(mu_filename),'" u 1:3 title "Im{mur}" w lp'

      if (.NOT.write_all_material_info_to_file) then
! only plot to the screen if we are operating in the interactive mode
        write(local_file_unit,'(A)')'set term wxt 0'
        write(local_file_unit,'(A,A,A)')'plot "',trim(mu_filename),'" u 1:2 title "Re{mur}" w lp,\'
        write(local_file_unit,'(A,A,A)')'     "',trim(mu_filename),'" u 1:3 title "Im{mur}" w lp'
        write(local_file_unit,'(A)')'pause 100'
      end if
      
      CLOSE(unit=local_file_unit)

! plot permeability response

      command='gnuplot '//trim(gnuplot_filename)//' & '
      CALL system(command)

    end if
  
  end do ! next material to output

  write(*,*)'FINISHED: Volume_material_frequency_response'
  
  RETURN

  
END SUBROUTINE Volume_material_frequency_response
