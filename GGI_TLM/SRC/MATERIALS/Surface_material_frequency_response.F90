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
! Name surface_material_frequency_response
!     
!
! Description:
!      Calculate the frequency response of surface materials
!      
!
! Comments:
!      
!      
!
! History
!
!     started 07/11/2012 CJS
!    3/09/2014		CJS: Implement simple diode impedance boundary conditions
!     17/10/2014 CJS - add write_all_material_info_to_file flag - eliminates the need for user input in this case
!

SUBROUTINE surface_material_frequency_response(write_all_material_info_to_file)

USE TLM_general
USE TLM_surface_materials
USE filter_types
USE filter_operators
USE filter_functions
USE file_information
USE file_general
USE constants

IMPLICIT NONE

! variables passed to subroutine
  
  logical :: write_all_material_info_to_file

! local_variables
  
  integer :: material_number,first_material,last_material
  character(LEN=filename_length)    :: z11_filename
  character(LEN=filename_length)    :: z12_filename
  character(LEN=filename_length)    :: z21_filename
  character(LEN=filename_length)    :: z22_filename
  real*8  :: fmin,fmax,fstep
  real*8  :: sigma
  character(LEN=filename_length)    :: gnuplot_filename
  character(LEN=filename_length)    :: base_filename
  character(LEN=filename_length)    :: jpeg_filename
  character(LEN=line_length)	    :: command

  integer	:: pol,pol1,pol2
  character*15     :: polstring(3)

! function_types

! START

  write(*,*)'CALLED: surface_material_frequency_response'
  write(*,*)
  write(*,*)'Number of surface materials=',n_surface_materials
  
  if (n_surface_materials.EQ.0) then
    write(*,*)'No surface materials to view'
    RETURN
  end if

! read the material number to output

! read the surface material number to view  
  if (write_all_material_info_to_file) then
! set the material number to zero to indicate that all material data should be written
    material_number=0
  else
    write(*,*)
    write(*,*)'Enter the surface material number to view or 0 to view all of them'
    read(*,*)material_number
  end if
  
  if ( (material_number.lt.0).OR.(material_number.gt.n_surface_materials) ) then
    write(*,*)'surface_material_number is outside the available range'
    write(*,*)'Number of surface materials=',n_surface_materials
    RETURN
  end if
  
  if (material_number.eq.0) then
    first_material=1
    last_material=n_surface_materials
  else
    first_material=material_number
    last_material=material_number
  end if

  do material_number=first_material,last_material

    if (surface_material_list(material_number)%type.EQ.surface_material_type_PEC) then
    
      write(*,*)'Material number',material_number,' is PEC'
      RETURN
  
    else if (surface_material_list(material_number)%type.EQ.surface_material_type_PMC) then
  
      write(*,*)'Material number',material_number,' is PMC'
      RETURN
  
    else if (surface_material_list(material_number)%type.EQ.surface_material_type_DIODE) then
  
      write(*,*)'Material number',material_number,' is DIODE'
      RETURN
  
    else if ( (surface_material_list(material_number)%type.EQ.surface_material_type_DISPERSIVE).OR. 	&
              (surface_material_list(material_number)%type.EQ.surface_material_type_ANISOTROPIC_DISPERSIVE) )then
    
      fmin=surface_material_list(material_number)%fmin
      fmax=surface_material_list(material_number)%fmax
      fstep=(fmax-fmin)/200d0
    
      sigma=0d0
    
      if (surface_material_list(material_number)%type.EQ.surface_material_type_DISPERSIVE) then
        pol1=1
        pol2=1
        polstring(1)=''
      else if (surface_material_list(material_number)%type.EQ.surface_material_type_ANISOTROPIC_DISPERSIVE) then
        pol1=1
        pol2=3   
        polstring(1)=':X_polarisation'
        polstring(2)=':Y_polarisation'
        polstring(3)=':Z_polarisation'
      end if 
    
      do pol=pol1,pol2

! Write Z11 filter response
  
        z11_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'.z11.fout'

! get the jpg output filenames - if we are running automatically then tag with the material number
        if (write_all_material_info_to_file) then
          base_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z11.jpg'
          CALL add_integer_to_filename(base_filename,material_number,jpeg_filename)
        else
          jpeg_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z11.jpg'
        end if
    
        CALL output_material_frequency_response(surface_material_list(material_number)%Z11_S(pol),	&
                                              sigma,fmin,fmax,fstep,local_file_unit,z11_filename)
      
! Write Z11 gnuplot file

        gnuplot_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z11.plt'
        OPEN(unit=local_file_unit,file=gnuplot_filename)
  
        write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
        write(local_file_unit,'(A)')'set ylabel "Z11 (ohms)"'
        write(local_file_unit,'(A,A,A,A)')'set title "',trim(surface_material_list(material_number)%name),trim(polstring(pol)),'"'

        write(local_file_unit,'(A)')'set term jpeg'
        write(local_file_unit,'(A,A,A)')'set output "',trim(jpeg_filename),'"'
        write(local_file_unit,'(A,A,A)')'plot "',trim(z11_filename),'" u 1:2 title "Re{z11}" w lp,\'
        write(local_file_unit,'(A,A,A)')'     "',trim(z11_filename),'" u 1:3 title "Im{z11}" w lp'

        if (.NOT.write_all_material_info_to_file) then
! only plot to the screen if we are operating in the interactive mode
          write(local_file_unit,'(A)')'set term wxt 0'
          write(local_file_unit,'(A,A,A)')'plot "',trim(z11_filename),'" u 1:2 title "Re{z11}" w lp,\'
          write(local_file_unit,'(A,A,A)')'     "',trim(z11_filename),'" u 1:3 title "Im{z11}" w lp'
          write(local_file_unit,'(A)')'pause 100'
        end if
	
        CLOSE(unit=local_file_unit)

! plot Z11 response

        command='gnuplot '//trim(gnuplot_filename)//' & '
        CALL system(command)

! Write Z12 filter response
  
        z12_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'.z12.fout'

! get the jpg output filenames - if we are running automatically then tag with the material number
        if (write_all_material_info_to_file) then
          base_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z12.jpg'
          CALL add_integer_to_filename(base_filename,material_number,jpeg_filename)
        else
          jpeg_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z12.jpg'
        end if
    
        CALL output_material_frequency_response(surface_material_list(material_number)%Z12_S(pol),	&
                                              sigma,fmin,fmax,fstep,local_file_unit,z12_filename)
! Write Z12 gnuplot file

        gnuplot_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z12.plt'
        OPEN(unit=local_file_unit,file=gnuplot_filename)
  
        write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
        write(local_file_unit,'(A)')'set ylabel "z12 (ohms)"'
        write(local_file_unit,'(A,A,A,A)')'set title "',trim(surface_material_list(material_number)%name),trim(polstring(pol)),'"'

        write(local_file_unit,'(A)')'set term jpeg'
        write(local_file_unit,'(A,A,A)')'set output "',trim(jpeg_filename),'"'
        write(local_file_unit,'(A,A,A)')'plot "',trim(z12_filename),'" u 1:2 title "Re{z12}" w lp,\'
        write(local_file_unit,'(A,A,A)')'     "',trim(z12_filename),'" u 1:3 title "Im{z12}" w lp'

        if (.NOT.write_all_material_info_to_file) then
! only plot to the screen if we are operating in the interactive mode
          write(local_file_unit,'(A)')'set term wxt 0'
          write(local_file_unit,'(A,A,A)')'plot "',trim(z12_filename),'" u 1:2 title "Re{z12}" w lp,\'
          write(local_file_unit,'(A,A,A)')'     "',trim(z12_filename),'" u 1:3 title "Im{z12}" w lp'
          write(local_file_unit,'(A)')'pause 100'
        end if
	
        CLOSE(unit=local_file_unit)

! plot Z12 response

        command='gnuplot '//trim(gnuplot_filename)//' & '
        CALL system(command)
      
! Write Z21 filter response  
  
        z21_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'.z21.fout'

! get the jpg output filenames - if we are running automatically then tag with the material number
        if (write_all_material_info_to_file) then
          base_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z21.jpg'
          CALL add_integer_to_filename(base_filename,material_number,jpeg_filename)
        else
          jpeg_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z21.jpg'
        end if
    
        CALL output_material_frequency_response(surface_material_list(material_number)%Z21_S(pol),	&
                                              sigma,fmin,fmax,fstep,local_file_unit,z21_filename)
! Write Z21 gnuplot file

        gnuplot_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z21.plt'
        OPEN(unit=local_file_unit,file=gnuplot_filename)
  
        write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
        write(local_file_unit,'(A)')'set ylabel "z21 (ohms)"'
        write(local_file_unit,'(A,A,A,A)')'set title "',trim(surface_material_list(material_number)%name),trim(polstring(pol)),'"'

        write(local_file_unit,'(A)')'set term jpeg'
        write(local_file_unit,'(A,A,A)')'set output "',trim(jpeg_filename),'"'
        write(local_file_unit,'(A,A,A)')'plot "',trim(z21_filename),'" u 1:2 title "Re{z21}" w lp,\'
        write(local_file_unit,'(A,A,A)')'     "',trim(z21_filename),'" u 1:3 title "Im{z21}" w lp'

        if (.NOT.write_all_material_info_to_file) then
! only plot to the screen if we are operating in the interactive mode
          write(local_file_unit,'(A)')'set term wxt 0'
          write(local_file_unit,'(A,A,A)')'plot "',trim(z21_filename),'" u 1:2 title "Re{z21}" w lp,\'
          write(local_file_unit,'(A,A,A)')'     "',trim(z21_filename),'" u 1:3 title "Im{z21}" w lp'
          write(local_file_unit,'(A)')'pause 100'
        end if
	
        CLOSE(unit=local_file_unit)

! plot Z21 response

        command='gnuplot '//trim(gnuplot_filename)//' & '
        CALL system(command)

! Write Z22 filter response 
  
        z22_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'.z22.fout'

! get the jpg output filenames - if we are running automatically then tag with the material number
        if (write_all_material_info_to_file) then
          base_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z22.jpg'
          CALL add_integer_to_filename(base_filename,material_number,jpeg_filename)
        else
          jpeg_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z22.jpg'
        end if
     
        CALL output_material_frequency_response(surface_material_list(material_number)%Z22_S(pol),	&
                                            sigma,fmin,fmax,fstep,local_file_unit,z22_filename)
! Write Z22 gnuplot file

        gnuplot_filename=trim(surface_material_list(material_number)%name)//trim(polstring(pol))//'_z22.plt'
        OPEN(unit=local_file_unit,file=gnuplot_filename)
  
        write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
        write(local_file_unit,'(A)')'set ylabel "z22 (ohms)"'
        write(local_file_unit,'(A,A,A,A)')'set title "',trim(surface_material_list(material_number)%name),trim(polstring(pol)),'"'

        write(local_file_unit,'(A)')'set term jpeg'
        write(local_file_unit,'(A,A,A)')'set output "',trim(jpeg_filename),'"'
        write(local_file_unit,'(A,A,A)')'plot "',trim(z22_filename),'" u 1:2 title "Re{z22}" w lp,\'
        write(local_file_unit,'(A,A,A)')'     "',trim(z22_filename),'" u 1:3 title "Im{z22}" w lp'

       if (.NOT.write_all_material_info_to_file) then
! only plot to the screen if we are operating in the interactive mode
          write(local_file_unit,'(A)')'set term wxt 0'
          write(local_file_unit,'(A,A,A)')'plot "',trim(z22_filename),'" u 1:2 title "Re{z22}" w lp,\'
          write(local_file_unit,'(A,A,A)')'     "',trim(z22_filename),'" u 1:3 title "Im{z22}" w lp'
          write(local_file_unit,'(A)')'pause 100'
        end if
	
        CLOSE(unit=local_file_unit)

! plot Z22 response

        command='gnuplot '//trim(gnuplot_filename)//' & '
        CALL system(command)
      
      end do ! next polarisation
    
    end if

  end do ! next material

  write(*,*)'FINISHED: surface_material_frequency_response'
  
  RETURN

  
END SUBROUTINE surface_material_frequency_response
