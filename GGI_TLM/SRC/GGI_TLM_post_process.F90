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
! PROGRAM GGI_TLM_post_process
!
! NAME
!     GGI_TLM_post_process
!
! DESCRIPTION
!     GGI_TLM post processor
!     Read and process TLM solver output files
!
!     Perform operations such as Fourier Transform, combination of results
!     Plot results with Gnuplot
!     Create animation files for visualisation with ParaView
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/08/2012 CJS
!
!
PROGRAM GGI_TLM_post_process

USE post_process
USE TLM_general
USE File_information

IMPLICIT NONE

! local variables

  integer,parameter	:: number_of_options=18
  integer	:: option

  character(len=256)	:: command

! START
  
  CALL write_progress('STARTED: GGI_TLM_post_process')

  CALL write_line('GGI_TLM_post_process',0,output_to_screen_flag)
  
  CALL write_license()
  
  OPEN(unit=record_user_inputs_unit,file='GGI_TLM_post_process_in_temp.txt')

10 CONTINUE  ! start of post_processing action

  write(*,*)
  write(*,*)'Post processing options are:'
  write(*,*)
  write(*,*)'1.  Extract Time Domain Data'
  write(*,*)'2.  Plot Time Domain Data'
  write(*,*)'3.  Fourier Transform (Time Domain to Frequency Domain)'
  write(*,*)'4.  Plot Frequency Domain Data'
  write(*,*)'5.  Create Animation of Time Domain Surface/Volume Field/ Current Output'
  write(*,*)'6.  Combine Frequency Domain Data'
  write(*,*)'7.  Combine Frequency Domain Magnitude Data'
  write(*,*)'8.  Frequency average'
  write(*,*)'9.  Oneport analysis'
  write(*,*)'10. Twoport analysis'
  write(*,*)'11. IELF analysis'
  write(*,*)'12. Create Animation of Frequency Domain Surface Field Output'
  write(*,*)'13. Extract mode from Frequency Domain Surface Field Output'
  write(*,*)'14. Sum Frequency Domain Data'
  write(*,*)'15. Transmission cross section calculation'
  write(*,*)'16. Fourier Transform with frequency warping (Time Domain to Frequency Domain)'
  write(*,*)'17. Create vector field Animation of Frequency Domain Surface Field Output'
  write(*,*)'18. Create vector field Animation of Frequency Domain Volume Field Output'
  write(*,*)
  
  write(*,'(A,I2,A)')'Please enter the required post processing option 1 :',number_of_options,' or 0 to quit'
  read(*,*)option
  
  if (option.EQ.0) then  ! close files, deallocate memory and stop
  
    write(record_user_inputs_unit,*)option,' Post processing option: Quit'
    CLOSE(unit=record_user_inputs_unit)
    
    command='mv GGI_TLM_post_process_in_temp.txt GGI_TLM_post_process_in.txt'
    CALL system(command)
  
    write(*,*)
    write(*,*)'The action of the GGI_TLM_post_processor just performed has been'
    write(*,*)"recorded in the file 'GGI_TLM_post_process_in.txt' "
    write(*,*)'The process can be re-run with the command:'
    write(*,*)'GGI_TLM_post_process < GGI_TLM_post_process_in.txt'
    write(*,*)
  
    CALL write_progress('FINISHED: GGI_TLM_post_process')
  
    STOP
    
  else if (option.EQ.1) then
  
    write(*,*)'Extract Time Domain Data'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: EXTRACT TIME DOMAIN DATA'
    CALL extract_Time_Domain_Data()
    
  else if (option.EQ.2) then
  
    write(*,*)'Plot Time Domain Data'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: PLOT TIME DOMAIN DATA'
    CALL plot_Time_Domain_Data()
    
  else if (option.EQ.3) then
    
    write(*,*)'Fourier Transform (Time Domain to Frequency Domain)'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: FOURIER TRANSFORM'
    CALL Fourier_Transform()

  else if (option.EQ.4) then
  
    write(*,*)'Plot Frequency Domain Data'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: PLOT FREQUENCY DOMAIN DATA'
    CALL plot_Frequency_Domain_Data()
    
  else if (option.EQ.5) then
  
    write(*,*)'Create Animation of Time Domain Surface/Volume Field/ Current Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE SURFACE/VOLUME ANIMATION'
    CALL create_animation()
   
  else if (option.EQ.6) then
  
    write(*,*)'Combine Frequency Domain Data'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: COMBINE FREQUENCY DOMAIN DATA'
    CALL Combine_Frequency_Domain_Data()
   
  else if (option.EQ.7) then
  
    write(*,*)'Combine Frequency Domain Magnitude Data'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: COMBINE FREQUENCY DOMAIN MAGNITUDE DATA'
    CALL Combine_Frequency_Domain_Magnitude_Data()
   
  else if (option.EQ.8) then
  
    write(*,*)'Frequency average'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: FREQUENCY AVERAGE'
    CALL Frequency_Average()
   
  else if (option.EQ.9) then
  
    write(*,*)'Oneport analysis'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: ONEPORT ANALYSIS'
    CALL Oneport()
   
  else if (option.EQ.10) then
  
    write(*,*)'Twoport analysis'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: TWOPORT ANALYSIS'
    CALL Twoport()
   
  else if (option.EQ.11) then
  
    write(*,*)'IELF analysis'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: IELF ANALYSIS'
    CALL IELF()
    
  else if (option.EQ.12) then
  
    write(*,*)'Create Animation of Frequency Domain Surface Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE FREQUENCY DOMAIN ANIMATION'
    CALL create_frequency_domain_animation()
    
  else if (option.EQ.13) then
  
    write(*,*)'Extract mode from Frequency Domain Surface Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: EXTRACT MODE FROM FREQUENCY DOMAIN OUTPUT'
    CALL extract_mode()
   
  else if (option.EQ.14) then
  
    write(*,*)'Sum Frequency Domain Data'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: SUM FREQUENCY DOMAIN DATA'
    CALL Sum_Frequency_Domain_Data()
   
  else if (option.EQ.15) then
  
    write(*,*)'Transmission Cross Section calculation'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: TRANSMISSION CROSS SECTION CALCULATION'
    CALL Transmission_Cross_Section()
    
  else if (option.EQ.16) then
    
    write(*,*)'Fourier Transform with frequency warping (Time Domain to Frequency Domain)'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: FOURIER TRANSFORM WITH FREQUENCY WARPING'
    CALL Fourier_Transform_warp()
    
  else if (option.EQ.17) then
  
    write(*,*)'Create Vector Field Animation of Frequency Domain Surface Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE VECTOR SURFACE FREQUENCY DOMAIN ANIMATION'
    CALL create_vector_surface_frequency_domain_animation()
    
  else if (option.EQ.18) then
  
    write(*,*)'Create Vector Field Animation of Frequency Domain Volume Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE VECTOR VOLUME FREQUENCY DOMAIN ANIMATION'
    CALL create_vector_volume_frequency_domain_animation()

  end if
  
  GOTO 10
  
END PROGRAM GGI_TLM_post_process
