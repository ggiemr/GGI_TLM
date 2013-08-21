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
! NAME
!     PROGRAM GGI_TLM_filter_format_convert
!
! DESCRIPTION
!     Given a material filter file in rational function form, read and convert the file
!     into rational function and pole-zero formats
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 8/1/2013 CJS
!
!
  PROGRAM GGI_TLM_filter_format_convert
!
       
USE TLM_general
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

  integer	:: function_loop
  
  real*8	:: fmin,fmax
  real*8	:: fstep,freq
  
  integer,parameter	:: n_frequency_values=100
  
  integer	:: freq_loop
  
  character*256	:: ip_line
  
  complex*16	:: fit_value_S
  complex*16	:: fit_value_S_PZ
  complex*16	:: fit_value_S_PR
  
  character(len=256) :: filter_filename
  character(len=256) :: filename(1:4)
 
! START

  CALL write_progress('STARTED: GGI_TLM_filter_format_convert')

  CALL write_line('GGI_TLM_filter_format_convert',0,output_to_screen_flag)
  
  CALL write_license()

! Read the filter name and model type 

  write(*,*)'Enter the name of the filter model file'
  read(*,'(A256)')FF_name
 
  CALL write_line('Filter model name:'//trim(FF_name),0,ff_output_to_screen)

  write(*,*)'Enter model type:'
  write(*,*)'1 : dielectric material model fit '
  write(*,*)'2 : magnetic material model fit '
  write(*,*)'3 : thin layer model fit'
  write(*,*)'4 : impedance model fit '
  
  read(*,*)fit_type
  
  if ( (fit_type.lt.1).or.(fit_type.gt.4) ) then
    write(*,*)'Error specifying fit_type'
    write(*,*)'Fit_type=',fit_type
    write(*,*)   &
    'Fittype must be 1 to 4'
    STOP
  end if

! Get material filter filename

  if (fit_type.eq.dielectric_material) then 
    filter_filename=trim(FF_name)//dielectric_filter_extension
    n_functions=2
  else if (fit_type.eq.magnetic_material) then  
    filter_filename=trim(FF_name)//magnetic_filter_extension
    n_functions=2
  else if (fit_type.eq.thin_layer) then  
    filter_filename=trim(FF_name)//thin_layer_filter_extension
    n_functions=4
  else if (fit_type.eq.impedance) then  
    filter_filename=trim(FF_name)//impedance_filter_extension
    n_functions=1
  end if
  
! allocate filters  
  ALLOCATE( filter_S(1:n_functions) )
  ALLOCATE( filter_S_PR(1:n_functions) )
  ALLOCATE( filter_S_PZ(1:n_functions) )
  ALLOCATE( filter_sigma(1:n_functions) )
  
! read filters 

  write(*,*)'Read filter(s)'

  OPEN(unit=filter_file_unit,file=filter_filename)
  
  read(filter_file_unit,'(A)')ip_line  ! filter title

  read(filter_file_unit,*)fmin,fmax
  
  if ( (fit_type.eq.dielectric_material).OR.	&
       (fit_type.eq.magnetic_material)       ) then 
    
    read(filter_file_unit,*)	! permittivity filter comment line
    CALL read_Sfilter(filter_S(1),filter_file_unit)
    read(filter_file_unit,*)filter_sigma(1)
    
    read(filter_file_unit,*)	! permeability filter comment line
    CALL read_Sfilter(filter_S(2),filter_file_unit)
    read(filter_file_unit,*)filter_sigma(2)
        
  else if (fit_type.eq.thin_layer) then  
    
    read(filter_file_unit,*)	! z11 filter comment line
    CALL read_Sfilter(filter_S(1),filter_file_unit)   
    read(filter_file_unit,*)	! z12 filter comment line
    CALL read_Sfilter(filter_S(2),filter_file_unit)    
    read(filter_file_unit,*)	! z21 filter comment line
    CALL read_Sfilter(filter_S(3),filter_file_unit)    
    read(filter_file_unit,*)	! z22 filter comment line
    CALL read_Sfilter(filter_S(4),filter_file_unit)
       
  else if (fit_type.eq.impedance) then  
    
    read(filter_file_unit,*)	! Impedance filter comment line
    CALL read_Sfilter(filter_S(1),filter_file_unit)
    
  end if
  
  CLOSE(unit=filter_file_unit)
 
! convert from rational function to pole-zero and pole residue formats

  write(*,*)'Convert from rational function to pole-zero and pole residue formats'

  do function_loop=1,n_functions
    filter_S_PZ(function_loop)=Convert_filter_S_to_S_PZ( filter_S(function_loop) )
    filter_S_PR(function_loop)=Convert_filter_S_to_S_PR( filter_S(function_loop) )
  end do ! next function

! write pole-zero filter function

  write(*,*)'Write pole-zero filter function'

  OPEN(unit=filter_file_unit,file=trim(filter_filename)//'_PZ')
  
  write(filter_file_unit,*)trim(ip_line)

  write(filter_file_unit,*)fmin,fmax
  
  if ( (fit_type.eq.dielectric_material).OR.	&
       (fit_type.eq.magnetic_material)       ) then 
    
    write(filter_file_unit,*)'permittivity filter'
    CALL write_S_PZ_filter(filter_S_PZ(1),filter_file_unit)
    write(filter_file_unit,*)filter_sigma(1),' electric conductivity'
    
    write(filter_file_unit,*)'permeability filter'
    CALL write_S_PZ_filter(filter_S_PZ(2),filter_file_unit)
    write(filter_file_unit,*)filter_sigma(2),' magnetic conductivity'
        
  else if (fit_type.eq.thin_layer) then  
    
    write(filter_file_unit,*)'z11 filter'
    CALL write_S_PZ_filter(filter_S_PZ(1),filter_file_unit)   
    write(filter_file_unit,*)'z12 filter'
    CALL write_S_PZ_filter(filter_S_PZ(2),filter_file_unit)    
    write(filter_file_unit,*)'z21 filter'
    CALL write_S_PZ_filter(filter_S_PZ(3),filter_file_unit)    
    write(filter_file_unit,*)'z22 filter'
    CALL write_S_PZ_filter(filter_S_PZ(4),filter_file_unit)
       
  else if (fit_type.eq.impedance) then  
    
    write(filter_file_unit,*)'Impedance filter'
    CALL write_S_PZ_filter(filter_S_PZ(1),filter_file_unit)
    
  end if
  
  CLOSE(unit=filter_file_unit)

! write pole-residue filter function

  write(*,*)'Write pole-residue filter function'

  OPEN(unit=filter_file_unit,file=trim(filter_filename)//'_PR')
  
  write(filter_file_unit,*)trim(ip_line)

  write(filter_file_unit,*)fmin,fmax
  
  if ( (fit_type.eq.dielectric_material).OR.	&
       (fit_type.eq.magnetic_material)       ) then 
    
    write(filter_file_unit,*)'permittivity filter'
    CALL write_S_PR_filter(filter_S_PR(1),filter_file_unit)
    write(filter_file_unit,*)filter_sigma(1),' electric conductivity'
    
    write(filter_file_unit,*)'permeability filter'
    CALL write_S_PR_filter(filter_S_PR(2),filter_file_unit)
    write(filter_file_unit,*)filter_sigma(2),' magnetic conductivity'
        
  else if (fit_type.eq.thin_layer) then  
    
    write(filter_file_unit,*)'z11 filter'
    CALL write_S_PR_filter(filter_S_PR(1),filter_file_unit)   
    write(filter_file_unit,*)'z12 filter'
    CALL write_S_PR_filter(filter_S_PR(2),filter_file_unit)    
    write(filter_file_unit,*)'z21 filter'
    CALL write_S_PR_filter(filter_S_PR(3),filter_file_unit)    
    write(filter_file_unit,*)'z22 filter'
    CALL write_S_PR_filter(filter_S_PR(4),filter_file_unit)
       
  else if (fit_type.eq.impedance) then  
    
    write(filter_file_unit,*)'Impedance filter'
    CALL write_S_PR_filter(filter_S_PR(1),filter_file_unit)
    
  end if
  
  CLOSE(unit=filter_file_unit)

! write filter function frequency responses

  write(*,*)'Write filter function frequency responses'

  if (fit_type.eq.dielectric_material) then 

    filename(1)=trim(filter_filename)//dielectric_trial_extension
    filename(2)=trim(filter_filename)//magnetic_trial_extension
 
  else if (fit_type.eq.magnetic_material) then 

    filename(1)=trim(filter_filename)//dielectric_trial_extension
    filename(2)=trim(filter_filename)//magnetic_trial_extension 
 
  else if (fit_type.eq.thin_layer) then 
  
    filename(1)=trim(filter_filename)//thin_layer_z11_trial_extension
    filename(2)=trim(filter_filename)//thin_layer_z12_trial_extension
    filename(3)=trim(filter_filename)//thin_layer_z21_trial_extension
    filename(4)=trim(filter_filename)//thin_layer_z22_trial_extension
 
  else if (fit_type.eq.impedance) then 
   
    filename(1)=trim(filter_filename)//impedance_trial_extension

  end if
  
  fstep=(fmax-fmin)/(n_frequency_values-1)
  
  do function_loop=1,n_functions
  
    OPEN(unit=filter_fit_data_file_unit,file=filename(function_loop))
  
    do freq_loop=1,n_frequency_values
    
      freq=fmin+(freq_loop-1)*fstep

! filter response    
      fit_value_S=evaluate_Sfilter_frequency_response(filter_S(function_loop),freq)
      fit_value_S_PZ=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(function_loop),freq)
      fit_value_S_PR=evaluate_Sfilter_PR_frequency_response(filter_S_PR(function_loop),freq)
      
! add the conductivity response if required
      if (fit_type.eq.dielectric_material) then 
  
        fit_value_S   =fit_value_S   +filter_sigma(1)/(j*2d0*pi*freq*eps0)
        fit_value_S_PZ=fit_value_S_PZ+filter_sigma(1)/(j*2d0*pi*freq*eps0)
        fit_value_S_PR=fit_value_S_PR+filter_sigma(1)/(j*2d0*pi*freq*eps0)
    
      end if 
        
      if (fit_type.eq.magnetic_material) then 
  
        fit_value_S   =fit_value_S   +filter_sigma(2)/(j*2d0*pi*freq*mu0)
        fit_value_S_PZ=fit_value_S_PZ+filter_sigma(2)/(j*2d0*pi*freq*mu0)
        fit_value_S_PR=fit_value_S_PR+filter_sigma(2)/(j*2d0*pi*freq*mu0)
    
      end if   
  
      write(filter_fit_data_file_unit,8000)	 freq	&
      						,real(fit_value_S)   ,imag(fit_value_S)	&
                                           	,real(fit_value_S_PZ),imag(fit_value_S_PZ)	&
					   	,real(fit_value_S_PR),imag(fit_value_S_PR)
8000  format(7E12.4)
     
    end do ! next freq_loop
    
  end do ! next function
  
! deallocate filters  
  DEALLOCATE( filter_S )
  DEALLOCATE( filter_S_PR )
  DEALLOCATE( filter_S_PZ )
  DEALLOCATE( filter_sigma )

  CALL write_progress('FINISHED: GGI_TLM_filter_format_convert')

  END
       
