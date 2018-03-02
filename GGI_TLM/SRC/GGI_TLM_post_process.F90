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

  integer,parameter	:: number_of_options=55
  integer	:: option

  character(len=256)	:: command

! START
  
  CALL write_progress('STARTED: GGI_TLM_post_process')

  CALL write_line('GGI_TLM_post_process',0,output_to_screen_flag)
  
  CALL write_license()
  
  OPEN(unit=record_user_inputs_unit,file='GGI_TLM_post_process_in_temp.txt')
  
  OPEN(unit=post_process_info_unit,file='GGI_TLM_post_process.info')
  write(post_process_info_unit,'(A)')"#START OF POST_PROCESSING DESCRIPTION"

10 CONTINUE  ! start of post_processing action

  write(*,*)
  write(*,*)'Post processing options are:'
  write(*,*)
  write(*,*)'1.  Extract Time Domain Data'
  write(*,*)'2.  Plot Time Domain Data'
  write(*,*)'3.  Discrete Time Fourier Transform (Time Domain to Frequency Domain)'
  write(*,*)'4.  Plot Frequency Domain Data'
  write(*,*)'5.  Create Animation of Time Domain Surface/Volume Field/ Current Output'
  write(*,*)'6.  Combine Frequency Domain Data: S(f)=A f1(f)/(B f2(f))+C'
  write(*,*)'7.  Combine Frequency Domain Magnitude Data: S(f)=A f1(f)/(B f2(f))+C'
  write(*,*)'8.  Frequency average'
  write(*,*)'9.  Oneport analysis'
  write(*,*)'10. Twoport analysis for symmetric structures'
  write(*,*)'11. IELF analysis'
  write(*,*)'12. Create Animation of Frequency Domain Surface Field Output'
  write(*,*)'13. Extract mode from Frequency Domain Surface Field Output'
  write(*,*)'14. Sum Frequency Domain Data'
  write(*,*)'15. Transmission cross section calculation'
  write(*,*)'16. Discrete Time Fourier Transform with frequency warping (Time Domain to Frequency Domain)'
  write(*,*)'17. Create vector field Animation of Frequency Domain Surface Field Output'
  write(*,*)'18. Create vector field Animation of Frequency Domain Volume Field Output'
  write(*,*)'19. S parameter to Z parameter transformation'
  write(*,*)'20. Calculate correlation matrix from time domain data'
  write(*,*)'21. Calculate correlation matrix from frequency domain data'
  write(*,*)'22. Calculate complex antenna factor from (S21) measurement for two identical antennas'
  write(*,*)'23. Calculate complex antenna factor from antenna coupling (S21) measurement with one unknown antenna'
  write(*,*)'24. Fast Fourier Transform (Time Domain to Frequency Domain)'
  write(*,*)'25. Generate random sample sequence either with or without an imposed correlation matrix'
  write(*,*)'26. Apply a filter function to time domain data'
  write(*,*)'27. Apply a filter function to frequency domain data'
  write(*,*)'28. Calculate the frequency domain transfer function of a filter function'
  write(*,*)'29. Combine Frequency Domain Data: S(f)=f1(f) f2(f)'
  write(*,*)'30. Calculate PDF and CDF from a data set'
  write(*,*)'31. Calculate correlation function from time domain data'
  write(*,*)'32. Calculate correlation function from frequency domain data'
  write(*,*)'33. Create Vector Animation of Time Domain Volume Field Output'
  write(*,*)'34. Create time domain near field scan'
  write(*,*)'35. Set random_number_seed'
  write(*,*)'36. Generate Far field plot data'
  write(*,*)'37. Re-format data file'
  write(*,*)'38. Visualise multi-dimensional data sets (up to 4D)'
  write(*,*)'39. Calculate the impulse response of a filter function'
  write(*,*)'40. S11 to VSWR'
  write(*,*)'41. Convert Frequency Domain Volume Field Output to x y z Re{field} Im{field} |field| format'
  write(*,*)'42. Calculate spatial correlation for 1D and 2D complex data sets'
  write(*,*)'43. Calculate and propagate Wigner function from complex spatial correlation data'
  write(*,*)'44. Sum Time Domain Data'
  write(*,*)'45. Square root Sum of squares/RMS calculation for Frequency Domain Data'
  write(*,*)'46. Reciprocal Frequency Domain Data'
  write(*,*)'47. Time-frequency analysis'
  write(*,*)'48. Statistical tools'
  write(*,*)'49. Create Poynting Vector or real/ imag E or H vector plot'
  write(*,*)'50. Multiply Time Domain Data'
  write(*,*)'51. Create Time Domain Force Vector Animation in conducting volumes'
  write(*,*)'52. Create_surface_frequency_domain_plot'
  write(*,*)'53. Calculate filter impulse response'
  write(*,*)'54. Create filter function (LPF, HPF)'
  write(*,*)'55. interpolate function(s)'
  write(*,*)
  
  write(*,'(A,I2,A)')'Please enter the required post processing option 1 :',number_of_options,' or 0 to quit'
  read(*,*)option
  
  if (option.EQ.0) then  ! close files, deallocate memory and stop
  
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: QUIT'
    CLOSE(unit=record_user_inputs_unit)
  
    write(post_process_info_unit,'(A)')"#END OF POST_PROCESSING DESCRIPTION"
    CLOSE(unit=post_process_info_unit)
    
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
    write(post_process_info_unit,*)'Extract time domain data'
    
    CALL extract_Time_Domain_Data()
    
  else if (option.EQ.2) then
  
    write(*,*)'Plot Time Domain Data'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: PLOT TIME DOMAIN DATA'
    write(post_process_info_unit,*)'Plot time domain data'
    CALL plot_Time_Domain_Data()
    
  else if (option.EQ.3) then
    
    write(*,*)'Fourier Transform (Time Domain to Frequency Domain)'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: FOURIER TRANSFORM'
    write(post_process_info_unit,*)'Fourier transform'
    CALL Fourier_Transform()

  else if (option.EQ.4) then
  
    write(*,*)'Plot Frequency Domain Data'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: PLOT FREQUENCY DOMAIN DATA'
    write(post_process_info_unit,*)'Plot frequency domain data'
    CALL plot_Frequency_Domain_Data()
    
  else if (option.EQ.5) then
  
    write(*,*)'Create Animation of Time Domain Surface/Volume Field/ Current Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE SURFACE/VOLUME ANIMATION'
    write(post_process_info_unit,*)'Create surface/volume animation'
    CALL create_animation()
   
  else if (option.EQ.6) then
  
    write(*,*)'Combine Frequency Domain Data'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: COMBINE FREQUENCY DOMAIN DATA: S(f)=A f1(f)/(B f2(f))+C'
    write(post_process_info_unit,*)'Combine frequency domain data: S(f)=A f1(f)/(B f2(f))+C'
    CALL Combine_Frequency_Domain_Data()
   
  else if (option.EQ.7) then
  
    write(*,*)'Combine Frequency Domain Magnitude Data'

    write(record_user_inputs_unit,*)option,	&
          ' POST PROCESSING OPTION: COMBINE FREQUENCY DOMAIN MAGNITUDE DATA: S(f)=A f1(f)/(B f2(f))+C'
    write(post_process_info_unit,*)'Combine frequency domain magnitude data: S(f)=A f1(f)/(B f2(f))+C'
    CALL Combine_Frequency_Domain_Magnitude_Data()
   
  else if (option.EQ.8) then
  
    write(*,*)'Frequency average'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: FREQUENCY AVERAGE'
    write(post_process_info_unit,*)'Frequency average'
    CALL Frequency_Average()
   
  else if (option.EQ.9) then
  
    write(*,*)'Oneport analysis'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: ONEPORT ANALYSIS'
    write(post_process_info_unit,*)'Oneport analysis'
    CALL Oneport()
   
  else if (option.EQ.10) then
  
    write(*,*)'Twoport analysis'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: TWOPORT ANALYSIS'
    write(post_process_info_unit,*)'Twoport analysis'
    CALL Twoport()
   
  else if (option.EQ.11) then
  
    write(*,*)'IELF analysis'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: IELF ANALYSIS'
    write(post_process_info_unit,*)'IELF analysis'
    CALL IELF()
    
  else if (option.EQ.12) then
  
    write(*,*)'Create Animation of Frequency Domain Surface Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE FREQUENCY DOMAIN ANIMATION'
    write(post_process_info_unit,*)'Create frequency domain animation'
    CALL create_frequency_domain_animation()
    
  else if (option.EQ.13) then
  
    write(*,*)'Extract mode from Frequency Domain Surface Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: EXTRACT MODE FROM FREQUENCY DOMAIN OUTPUT'
    write(post_process_info_unit,*)'Extract mode from frequency domain output'
    CALL extract_mode()
   
  else if (option.EQ.14) then
  
    write(*,*)'Sum Frequency Domain Data'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: SUM FREQUENCY DOMAIN DATA'
    write(post_process_info_unit,*)'Sum frequency domain data'
    CALL Sum_Frequency_Domain_Data()
   
  else if (option.EQ.15) then
  
    write(*,*)'Transmission Cross Section calculation'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: TRANSMISSION CROSS SECTION CALCULATION'
    write(post_process_info_unit,*)'Transmission cross section calculation'
    CALL Transmission_Cross_Section()
    
  else if (option.EQ.16) then
    
    write(*,*)'Fourier Transform with frequency warping (Time Domain to Frequency Domain)'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: FOURIER TRANSFORM WITH FREQUENCY WARPING'
    write(post_process_info_unit,*)'Fourier transform with frequency warping'
    CALL Fourier_Transform_warp()
    
  else if (option.EQ.17) then
  
    write(*,*)'Create Vector Field Animation of Frequency Domain Surface Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE VECTOR SURFACE FREQUENCY DOMAIN ANIMATION'
    write(post_process_info_unit,*)'Create vector surface frequency domain animation'
    CALL create_vector_surface_frequency_domain_animation()
    
  else if (option.EQ.18) then
  
    write(*,*)'Create Vector Field Animation of Frequency Domain Volume Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE VECTOR VOLUME FREQUENCY DOMAIN ANIMATION'
    write(post_process_info_unit,*)'Create vector volume frequency domain animation'
    CALL create_vector_volume_frequency_domain_animation()
    
  else if (option.EQ.19) then
  
    write(*,*)'S parameter to Z parameter transformation for symmetric structures'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: S PARAMETER TO Z PARAMETER TRANSFORMATION '
    write(post_process_info_unit,*)'S parameter to Z parameter transformation'
    CALL S_to_Z()
    
  else if (option.EQ.20) then
  
    write(*,*)'Calculate correlation of time domain data'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE CORRELATION MATRIX OF TIME DOMAIN DATA'
    write(post_process_info_unit,*)'Calculate correlation matrix of time domain data'
    CALL correlation_time()
    
  else if (option.EQ.21) then
  
    write(*,*)'Calculate correlation of frequency domain data'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE CORRELATION MATRIX OF FREQUENCY DOMAIN DATA'
    write(post_process_info_unit,*)'Calculate correlation matrix of frequency domain data'
    CALL correlation_frequency()
    
  else if (option.EQ.22) then
  
    write(*,*)'Calculate complex antenna factor from two antenna coupling (S21) measurement'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE COMPLEX ANTENNA FACTOR, TWO UNKNOWN ANTENNAS'
    write(post_process_info_unit,*)'Calculate complex antenna factor, two unknown antennas'
    CALL complex_antenna_factor()
    
  else if (option.EQ.23) then
  
    write(*,*)'Calculate complex antenna factor from two antenna coupling (S21) measurement'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE COMPLEX ANTENNA FACTOR, ONE UNKNOWN ANTENNA'
    write(post_process_info_unit,*)'Calculate complex antenna factor, one unknown antenna'
    CALL complex_antenna_factor_2()
    
  else if (option.EQ.24) then
  
    write(*,*)'Fast Fourier Transform (Time Domain to Frequency Domain)'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: FAST FOURIER TRANSFORM (TIME DOMAIN TO FREQUENCY DOMAIN)'
    write(post_process_info_unit,*)'Fast fourier transform (time domain to frequency domain)'
    CALL FFT_MAIN()
    
  else if (option.EQ.25) then
  
    write(*,*)'Generate random sample sequence either with or without an imposed correlation matrix'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: GENERATE RANDOM SAMPLE SEQUENCE'
    write(post_process_info_unit,*)'Generate random sample sequence'
    CALL Generate_random_sequence()
    
  else if (option.EQ.26) then
  
    write(*,*)'Apply a filter function to time domain data'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: APPLY A FILTER FUNCTION TO TIME DOMAIN DATA'
    write(post_process_info_unit,*)'Apply a filter function to time domain data'
    CALL Apply_filter_in_time_domain()
    
  else if (option.EQ.27) then
  
    write(*,*)'Apply a filter function to frequency domain data'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: APPLY A FILTER FUNCTION TO FREQUENCY DOMAIN DATA'
    write(post_process_info_unit,*)'Apply a filter function to frequency domain data'
    CALL Apply_filter_in_frequency_domain()
    
  else if (option.EQ.28) then
  
    write(*,*)'Calculate the frequency domain transfer function of a filter function'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE THE FREQUENCY DOMAIN TRANSFER FUNCTION OF A FILTER FUNCTION'
    write(post_process_info_unit,*)'Calculate the frequency domain transfer function of a filter function'
    CALL Calculate_filter_frequency_response()
    
  else if (option.EQ.29) then
  
    write(*,*)'Combine Frequency Domain Data: S(f)=f1(f) f2(f)'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: COMBINE FREQUENCY DOMAIN DATA: S(f)=F1(f) F2(f)'
    write(post_process_info_unit,*)'Combine frequency domain data: S(f)=F1(f) F2(f)'
    CALL Combine_Frequency_Domain_Data_2()
    
  else if (option.EQ.30) then
  
    write(*,*)'Calculate PDF and CDF from a data set'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE PDF AND CDF FROM A DATA SET'
    write(post_process_info_unit,*)'Calculate pdf and cdf from a data set'
    CALL PDF_CDF()
    
  else if (option.EQ.31) then
  
    write(*,*)'Calculate correlation of time domain data'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE CORRELATION FUNCTION OF TIME DOMAIN DATA'
    write(post_process_info_unit,*)'Calculate correlation function of time domain data'
    CALL correlation_function_time()
      
  else if (option.EQ.32) then
  
    write(*,*)'Calculate correlation of frequency domain data'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE CORRELATION FUNCTION OF FREQUENCY DOMAIN DATA'
    write(post_process_info_unit,*)'Calculate correlation function of frequency domain data'
    CALL correlation_function_frequency()
    
  else if (option.EQ.33) then
  
    write(*,*)'Create Vector Animation of Time Domain Volume Field Output'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE TIME DOMAIN VECTOR VOLUME ANIMATION'
    write(post_process_info_unit,*)'Create time domain vector volume animation'
    CALL create_time_domain_vector_animation()
    
  else if (option.EQ.34) then
  
    write(*,*)'Create time domain near field scan'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE_TIME_DOMAIN_NEAR_FIELD_SCAN'
    write(post_process_info_unit,*)'Create_time_domain_near_field_scan'
    CALL create_time_domain_near_field_scan()
    
  else if (option.EQ.35) then
  
    write(*,*)'Set random number seed'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: SET RANDOM NUMBER SEED'
    write(post_process_info_unit,*)'Set random number seed'
    CALL set_random_seed()
    
  else if (option.EQ.36) then
  
    write(*,*)'Generate Far field plot data'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: GENERATE FAR FIELD PLOT DATA'
    write(post_process_info_unit,*)'Generate far field plot data'
    CALL generate_far_field_plot()
    
  else if (option.EQ.37) then
  
    write(*,*)'Re-format data file'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: RE-FORMAT DATA FILE'
    write(post_process_info_unit,*)'Re-format data file'
    CALL re_format_data_file()
    
  else if (option.EQ.38) then
  
    write(*,*)'Visualise multi-dimensional data'
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: VISUALISE MULTI-DIMENSIONAL DATA'
    write(post_process_info_unit,*)'Visualise multi-dimensional data'
    CALL Vis_nD()
    
  else if (option.EQ.39) then
  
    write(*,*)'Calculate the impulse response of a filter function'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: CALCULATE THE IMPULSE RESPONSE OF A FILTER FUNCTION'
    write(post_process_info_unit,*)'Calculate the impulse response of a filter function'
    CALL Calculate_filter_time_response()
    
  else if (option.EQ.40) then
  
    write(*,*)'40. S11 to VSWR'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: S11 to VSWR '
    write(post_process_info_unit,*)'Convert S11 to VSWR'
    CALL S11_TO_VSWR()
    
  else if (option.EQ.41) then
  
    write(*,*)'41. Convert Frequency Domain Volume Field Output to x y z Re{field} Im{field} format'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: Convert Frequency Domain Volume Field Output '
    write(post_process_info_unit,*)'Convert Frequency Domain Volume Field Output to x y z Re{field} Im{field} format'
    CALL Convert_Frequency_Domain_Volume_Field_Format()
    
  else if (option.EQ.42) then
  
    write(*,*)'42. Calculate spatial correlation for 1D and 2D complex data sets'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: Calculate spatial correlation for 1D and 2D complex data sets '
    write(post_process_info_unit,*)'Calculate spatial correlation for 1D and 2D complex data sets'
    CALL Spatial_correlation()
    
  else if (option.EQ.43) then
  
    write(*,*)'43. Calculate and propagate Wigner function from complex spatial correlation data'
    write(record_user_inputs_unit,*)option,	&
    ' POST PROCESSING OPTION: Calculate and propagate Wigner function from complex spatial correlation data '
    write(post_process_info_unit,*)'Calculate and propagate Wigner function from complex spatial correlation data'
    CALL Wigner_function()
    
  else if (option.EQ.44) then
  
    write(*,*)'Sum Time Domain Data'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: SUM TIME DOMAIN DATA'
    write(post_process_info_unit,*)'Sum time domain data'
    CALL Sum_Time_Domain_Data()
   
  else if (option.EQ.45) then

    write(*,*)'Square root Sum of squares/RMS calculation for Frequency Domain Data'
    
    write(record_user_inputs_unit,*)option,	&
      ' POST PROCESSING OPTION: SQUARE ROOT SUM OF SQUARES SQUARES/RMS CALCULATION FOR FREQUENCY DOMAIN DATA'
    write(post_process_info_unit,*)'Square root Sum of squares/RMS calculation for Frequency Domain Data'
    CALL RMS_Frequency_Domain_Data()
    
  else if (option.EQ.46) then

    write(*,*)'Reciprocal Frequency Domain Data'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: RECIPROCAL FREQUENCY DOMAIN DATA'
    write(post_process_info_unit,*)'Reciprocal Frequency Domain Data'
    CALL Reciprocal_Frequency_Domain_Data()
    
  else if (option.EQ.47) then

    write(*,*)'Time-frequency analysis'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: Time-frequency analysis'
    write(post_process_info_unit,*)'Time-frequency analysis'
    CALL Time_frequency_analysis()
    
  else if (option.EQ.48) then

    write(*,*)'Statistical tools'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: STATISTICAL TOOLS'
    write(post_process_info_unit,*)'Apply statistical tools'
    CALL Statistical_tools()
    
  else if (option.EQ.49) then

    write(*,*)'Create Poynting Vector plot or real/ imag E or H vector plot'
    
    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE POYNTING VECTOR OR REAL/ IMAG E OR H VECTOR PLOT'
    write(post_process_info_unit,*)'Create Poynting Vector or real/ imag E or H vector plot'
    CALL create_poynting_vector_plot()
    
  else if (option.EQ.50) then
  
    write(*,*)'Multiply Time Domain Data'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: MULTIPLY TIME DOMAIN DATA'
    write(post_process_info_unit,*)'Multiply time domain data'
    CALL Multiply_Time_Domain_Data()
    
  else if (option.EQ.51) then
  
    write(*,*)'Create Time Domain Force Vector Animation in conducting volumes'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE TIME DOMAIN FORCE VECTOR ANIMATION'
    write(post_process_info_unit,*)'Create Time Domain Force Vector Animation in conducting volumes'
    CALL create_time_domain_force_vector_animation()
    
  else if (option.EQ.52) then
  
    write(*,*)'Create surface frequency domain plot'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE SURFACE FREQUENCY DOMAIN PLOT'
    write(post_process_info_unit,*)'Create surface frequency domain plot'
    CALL create_surface_frequency_domain_plot()
    
  else if (option.EQ.53) then
  
    write(*,*)'Calculate filter impulse response'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CALCULATE FILTER IMPULSE RESPONSE'
    write(post_process_info_unit,*)'Calculate filter impulse response'
    CALL filter_impulse_response()
    
  else if (option.EQ.54) then
  
    write(*,*)'Create filter function (LPF, HPF)'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: CREATE FILTER FUNCTION (LPF, HPF)'
    write(post_process_info_unit,*)'Create filter function (LPF, HPF)'
    CALL create_filter_function()
    
  else if (option.EQ.55) then
  
    write(*,*)'Interpolate function(s)'

    write(record_user_inputs_unit,*)option,' POST PROCESSING OPTION: INTERPOLATE FUNCTION(S)'
    write(post_process_info_unit,*)'Interpolate function(s)'
    CALL interpolate()
      
  else
   
    write(*,*)'Unknown option',option
    write(*,'(A,I2)')'Option should be in the range 0 to',number_of_options
    
    CALL write_progress('FAILED: GGI_TLM_post_process')
  
    STOP

  end if
  
  GOTO 10
  
END PROGRAM GGI_TLM_post_process
