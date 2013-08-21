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
!

SUBROUTINE surface_material_frequency_response()

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

! local_variables
  
  integer :: material_number
  character(LEN=filename_length)    :: z11_filename
  character(LEN=filename_length)    :: z12_filename
  character(LEN=filename_length)    :: z21_filename
  character(LEN=filename_length)    :: z22_filename
  real*8  :: fmin,fmax,fstep
  real*8  :: sigma
  character(LEN=filename_length)    :: gnuplot_filename
  character(LEN=line_length)	:: command

! function_types

! START

  write(*,*)'CALLED: surface_material_frequency_response'

! read the material number to output

! read the surface material number to view  
  write(*,*)
  write(*,*)'Enter the surface material number to view'
  read(*,*)material_number
  
  if ( (material_number.le.0).OR.(material_number.gt.n_surface_materials) ) then
    write(*,*)'surface_material_number is outside the available range'
    write(*,*)'Number of surface materials=',n_surface_materials
    RETURN
  end if

  if (surface_material_list(material_number)%type.EQ.surface_material_type_PEC) then
  
    write(*,*)'Material number',material_number,' is PEC'
    RETURN
  
  else if (surface_material_list(material_number)%type.EQ.surface_material_type_PMC) then
  
    write(*,*)'Material number',material_number,' is PMC'
    RETURN
  
  else if (surface_material_list(material_number)%type.EQ.surface_material_type_DISPERSIVE) then
    
    fmin=surface_material_list(material_number)%fmin
    fmax=surface_material_list(material_number)%fmax
    fstep=(fmax-fmin)/200d0
    
    sigma=0d0

! Write Z11 filter response
  
    z11_filename=trim(surface_material_list(material_number)%name)//'.z11.fout'
    
    CALL output_material_frequency_response(surface_material_list(material_number)%Z11_S,	&
                                            sigma,fmin,fmax,fstep,local_file_unit,z11_filename)
      
! Write Z11 gnuplot file

    gnuplot_filename=trim(surface_material_list(material_number)%name)//'_z11.plt'
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "Z11 (ohms)"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(surface_material_list(material_number)%name),'"'

    write(local_file_unit,'(A)')'set term jpeg'
    write(local_file_unit,'(A,A,A)')'set output "',trim(surface_material_list(material_number)%name)//'_z11.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z11_filename),'" u 1:2 title "Re{z11}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z11_filename),'" u 1:3 title "Im{z11}" w lp'

    write(local_file_unit,'(A)')'set term wxt 0'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z11_filename),'" u 1:2 title "Re{z11}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z11_filename),'" u 1:3 title "Im{z11}" w lp'
    write(local_file_unit,'(A)')'pause 100'
  
    CLOSE(unit=local_file_unit)

! plot Z11 response

    command='gnuplot '//trim(gnuplot_filename)//' & '
    CALL system(command)

! Write Z12 filter response
  
    z12_filename=trim(surface_material_list(material_number)%name)//'.z12.fout'
    
    CALL output_material_frequency_response(surface_material_list(material_number)%Z12_S,	&
                                            sigma,fmin,fmax,fstep,local_file_unit,z12_filename)
! Write Z12 gnuplot file

    gnuplot_filename=trim(surface_material_list(material_number)%name)//'_z12.plt'
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "z12 (ohms)"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(surface_material_list(material_number)%name),'"'

    write(local_file_unit,'(A)')'set term jpeg'
    write(local_file_unit,'(A,A,A)')'set output "',trim(surface_material_list(material_number)%name)//'_z12.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z12_filename),'" u 1:2 title "Re{z12}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z12_filename),'" u 1:3 title "Im{z12}" w lp'

    write(local_file_unit,'(A)')'set term wxt 0'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z12_filename),'" u 1:2 title "Re{z12}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z12_filename),'" u 1:3 title "Im{z12}" w lp'
    write(local_file_unit,'(A)')'pause 100'
  
    CLOSE(unit=local_file_unit)

! plot Z12 response

    command='gnuplot '//trim(gnuplot_filename)//' & '
    CALL system(command)
      
! Write Z21 filter response  
  
    z21_filename=trim(surface_material_list(material_number)%name)//'.z21.fout'
    
    CALL output_material_frequency_response(surface_material_list(material_number)%Z21_S,	&
                                            sigma,fmin,fmax,fstep,local_file_unit,z21_filename)
! Write Z21 gnuplot file

    gnuplot_filename=trim(surface_material_list(material_number)%name)//'_z21.plt'
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "z21 (ohms)"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(surface_material_list(material_number)%name),'"'

    write(local_file_unit,'(A)')'set term jpeg'
    write(local_file_unit,'(A,A,A)')'set output "',trim(surface_material_list(material_number)%name)//'_z21.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z21_filename),'" u 1:2 title "Re{z21}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z21_filename),'" u 1:3 title "Im{z21}" w lp'

    write(local_file_unit,'(A)')'set term wxt 0'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z21_filename),'" u 1:2 title "Re{z21}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z21_filename),'" u 1:3 title "Im{z21}" w lp'
    write(local_file_unit,'(A)')'pause 100'
  
    CLOSE(unit=local_file_unit)

! plot Z21 response

    command='gnuplot '//trim(gnuplot_filename)//' & '
    CALL system(command)

! Write Z22 filter response 
  
    z22_filename=trim(surface_material_list(material_number)%name)//'.z22.fout'
    
    CALL output_material_frequency_response(surface_material_list(material_number)%Z22_S,	&
                                            sigma,fmin,fmax,fstep,local_file_unit,z22_filename)
! Write Z22 gnuplot file

    gnuplot_filename=trim(surface_material_list(material_number)%name)//'_z22.plt'
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "z22 (ohms)"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(surface_material_list(material_number)%name),'"'

    write(local_file_unit,'(A)')'set term jpeg'
    write(local_file_unit,'(A,A,A)')'set output "',trim(surface_material_list(material_number)%name)//'_z22.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z22_filename),'" u 1:2 title "Re{z22}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z22_filename),'" u 1:3 title "Im{z22}" w lp'

    write(local_file_unit,'(A)')'set term wxt 0'
    write(local_file_unit,'(A,A,A)')'plot "',trim(z22_filename),'" u 1:2 title "Re{z22}" w lp,\'
    write(local_file_unit,'(A,A,A)')'     "',trim(z22_filename),'" u 1:3 title "Im{z22}" w lp'
    write(local_file_unit,'(A)')'pause 100'
  
    CLOSE(unit=local_file_unit)

! plot Z22 response

    command='gnuplot '//trim(gnuplot_filename)//' & '
    CALL system(command)

  end if

  write(*,*)'FINISHED: surface_material_frequency_response'
  
  RETURN

  
END SUBROUTINE surface_material_frequency_response
