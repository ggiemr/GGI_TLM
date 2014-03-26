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
! SUBROUTINE write_filter_fit_gnuplot_file
!
! NAME
!     write_filter_fit_gnuplot_file
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/12/2013 CJS
!
!
SUBROUTINE write_filter_fit_gnuplot_file

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff

IMPLICIT NONE
 
! local variables
  
  character(len=256) :: input_data_filename
  character(len=256) :: model_order_filename
  character(len=256) :: temp_filename
  character(len=256) :: title_filename
  character(len=256) :: model_fit_filename
  character(len=256) :: gnuplot_filename
  character(len=11)  :: stability_string
  
  character(len=256) :: title
  character(len=16)  :: error_string

! START

  CALL write_line('CALLED: write_filter_fit_gnuplot_file',0,ff_output_to_screen)
  
  CALL add_integer_to_filename(FF_name,order,model_order_filename)
  
  temp_filename=trim(FF_name)//' Model order '
  CALL add_integer_to_filename(temp_filename,order,title_filename)
  
  if (stable_filter) then
    stability_string=' : STABLE'
  else
    stability_string=' : UNSTABLE'
  end if 
  
  write(error_string,'(F16.6)')Mean_square_error

  if (fit_type.eq.dielectric_material) then 

! PLOT FILE FOR RELATIVE PERMITTIVITY DATA  

    input_data_filename=trim(FF_name)//dielectric_input_data_extension  
    model_fit_filename=trim(model_order_filename)//dielectric_trial_extension

    gnuplot_filename=trim(model_order_filename)//'_eps.plt'
    
    title=trim(title_filename)//trim(stability_string)//' : Error='//error_string
    
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "Relative Permittivity"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(title),'"'
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'_eps.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{epsr}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{epsr}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{epsr}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{epsr}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
  
    CLOSE(unit=local_file_unit)
 
  else if (fit_type.eq.magnetic_material) then 
  
    input_data_filename=trim(FF_name)//magnetic_input_data_extension  
    model_fit_filename=trim(model_order_filename)//magnetic_trial_extension

    gnuplot_filename=trim(model_order_filename)//'_mu.plt'
    
    title=trim(title_filename)//trim(stability_string)//' : Error='//error_string
    
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "Relative Permeability"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(title),'"'
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'_eps.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{mur}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{mur}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{mur}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{mur}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
  
    CLOSE(unit=local_file_unit)
 
  else if (fit_type.eq.thin_layer) then 

    gnuplot_filename=trim(model_order_filename)//'_Zij.plt'
    
    title=trim(title_filename)//trim(stability_string)//' : Error='//error_string
    
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    input_data_filename=trim(FF_name)//thin_layer_z11_input_data_extension
    model_fit_filename =trim(model_order_filename)//thin_layer_z11_trial_extension
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "Impedance (ohms)"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(title),'"'
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'_z11.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{z11}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{z11}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{z11}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{z11}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
    
    input_data_filename=trim(FF_name)//thin_layer_z12_input_data_extension
    model_fit_filename =trim(model_order_filename)//thin_layer_z12_trial_extension    
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'_z12.jpeg','"'
    write(local_file_unit,'(A)')'set autoscale x'
    write(local_file_unit,'(A)')'set autoscale y'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{z12}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{z12}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{z12}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{z12}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
    
    input_data_filename=trim(FF_name)//thin_layer_z21_input_data_extension
    model_fit_filename =trim(model_order_filename)//thin_layer_z21_trial_extension    
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'_z21.jpeg','"'
    write(local_file_unit,'(A)')'set autoscale x'
    write(local_file_unit,'(A)')'set autoscale y'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{z21}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{z21}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{z21}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{z21}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
    
    input_data_filename=trim(FF_name)//thin_layer_z22_input_data_extension
    model_fit_filename =trim(model_order_filename)//thin_layer_z22_trial_extension     
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'_z22.jpeg','"'
    write(local_file_unit,'(A)')'set autoscale x'
    write(local_file_unit,'(A)')'set autoscale y'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{z22}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{z22}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{z22}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{z22}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
  
    CLOSE(unit=local_file_unit)
 
  else if (fit_type.eq.impedance) then 
   
    input_data_filename=trim(FF_name)//impedance_input_data_extension
    model_fit_filename =trim(model_order_filename)//impedance_trial_extension

    gnuplot_filename=trim(model_order_filename)//'_Z.plt'
    
    title=trim(title_filename)//trim(stability_string)//' : Error='//error_string
    
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "Impedance (ohms)"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(title),'"'
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'_Z.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{Z}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{Z}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{Z}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{Z}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
  
    CLOSE(unit=local_file_unit)
 
  else if (fit_type.eq.general) then 
   
    input_data_filename=trim(FF_name)//general_input_data_extension
    model_fit_filename =trim(model_order_filename)//general_trial_extension

    gnuplot_filename=trim(model_order_filename)//'.plt'
    
    title=trim(title_filename)//trim(stability_string)//' : Error='//error_string
    
    OPEN(unit=local_file_unit,file=gnuplot_filename)
  
    write(local_file_unit,'(A)')'set xlabel "Frequency (Hz)"'
    write(local_file_unit,'(A)')'set ylabel "Function value"'
    write(local_file_unit,'(A,A,A)')'set title "',trim(title),'"'
    write(local_file_unit,'(A)')'#set term jpeg'
    write(local_file_unit,'(A,A,A)')'#set output "',trim(model_order_filename)//'.jpeg','"'
    write(local_file_unit,'(A,A,A)')'plot "',trim(input_data_filename),'" u 1:2 title "Re{F}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(input_data_filename),'" u 1:3 title "Im{F}: Input data" w p, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:2 title "Re{F}: Model fit " w l, \'
    write(local_file_unit,'(A,A,A)')'     "',trim(model_fit_filename) ,'" u 1:3 title "Im{F}: Model fit " w l'
    write(local_file_unit,'(A)')'pause -1'
  
    CLOSE(unit=local_file_unit)

  end if

  CALL write_line('FINISHED: write_filter_fit_gnuplot_file',0,ff_output_to_screen)

END SUBROUTINE write_filter_fit_gnuplot_file
