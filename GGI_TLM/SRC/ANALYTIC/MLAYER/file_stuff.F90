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
!  SUBROUTINE open_files
!  SUBROUTINE read_Mlayer_input_file
!  SUBROUTINE close_files

! NAME
!     SUBROUTINE  open_files
!
! DESCRIPTION
! 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
  SUBROUTINE open_files()
  
USE Mlayer_file_module
  
! START
  
  write(*,*)'Enter the input filename (without .inp extension)'
  read(*,'(A256)')problem_name
    
  write(*,*)'Opening file:',trim(problem_name)//input_file_extn

  open(UNIT=input_file_unit,						&
       FILE=trim(problem_name)//input_file_extn,	&
       STATUS='old',							&
       ERR=9000)
       
  open(UNIT=info_file_unit,						&
       FILE=trim(problem_name)//info_file_extn,	&
       ERR=9010)
       
  write(info_file_unit,*)'Problem name=',problem_name (1:problem_name_length) 

  write(info_file_unit,*)'Opening S parameter files'

  open(UNIT=S11_file_unit,						&
       FILE=trim(problem_name)//S11_file_extn,	&
       ERR=9020)
       
  open(UNIT=S12_file_unit,						&
       FILE=trim(problem_name)//S12_file_extn,	&
       ERR=9020)
       
  open(UNIT=S21_file_unit,						&
       FILE=trim(problem_name)//S21_file_extn,	&
       ERR=9020)
       
  open(UNIT=S22_file_unit,						&
       FILE=trim(problem_name)//S22_file_extn,	&
       ERR=9020)
       
  write(info_file_unit,*)'Opening Z parameter files'
  
  open(UNIT=Z11_file_unit,						&
       FILE=trim(problem_name)//Z11_file_extn,	&
       ERR=9030)
       
  open(UNIT=Z12_file_unit,						&
       FILE=trim(problem_name)//Z12_file_extn,	&
       ERR=9030)
       
  open(UNIT=Z21_file_unit,						&
       FILE=trim(problem_name)//Z21_file_extn,	&
       ERR=9030)
       
  open(UNIT=Z22_file_unit,						&
       FILE=trim(problem_name)//Z22_file_extn,	&
       ERR=9030)
       
  open(UNIT=power_file_unit,						&
       FILE=trim(problem_name)//power_file_extn,	&
       ERR=9030)

  return

9000 continue
  write(*,*)'Error opening input file'
  write(*,*)'Filename='
  write(*,*)trim(problem_name)//input_file_extn
  stop
  
9010 continue
  write(*,*)'Error opening info file'
  write(*,*)'Filename='
  write(*,*)trim(problem_name)//input_file_extn
  stop
  
9020 continue
  write(*,*)'Error opening S parameter file'
  stop
  
9030 continue
  write(*,*)'Error opening Z parameter file'
  stop

  END SUBROUTINE open_files

! NAME
!     SUBROUTINE  read_Mlayer_input_file
!
! DESCRIPTION
! 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
  SUBROUTINE read_Mlayer_input_file()
  
USE constants
USE Mlayer_module
USE Mlayer_file_module
USE filter_types
USE filter_functions

! local variables

  integer layer,layer_number
  character layer_type_in
  integer   material_name_length
  integer   thin_layer_name_length
  
  type(SFilter) :: filter_in
  
! START

  write(info_file_unit,*)' '
  write(info_file_unit,*)'Reading input file...'
  write(info_file_unit,*)' '
  
  read(input_file_unit,*,END=9000,ERR=9010)m
  write(*,*)'Number of layers=',m
  write(info_file_unit,*)'Number of layers=',m
  
  read(input_file_unit,*,END=9000,ERR=9010)angle
  write(*,*)'Incident angle=',angle
  write(info_file_unit,*)'Incident angle=',angle
  
  angle_rad=angle*pi/180d0
    
  read(input_file_unit,'(A2)',END=9000,ERR=9010)polarisation_in
  write(*,*)'Polarisation=',polarisation_in
  write(info_file_unit,*)'Polarisation=',polarisation_in
  
  polarisation=0
  if ( (polarisation_in(1:2).EQ."TE").OR.(polarisation_in(1:2).EQ."te") ) then
    polarisation=te
  else if ( (polarisation_in(1:2).EQ."TM").OR.(polarisation_in(1:2).EQ."tm") ) then
    polarisation=tm
  else
    write(*,*)'Unknown polarisation:',polarisation_in
    write(*,*)'The polarisation should be either TE or TM'
  end if
  
  write(*,*)'Polarisation type=',polarisation
  write(info_file_unit,*)'Polarisation type=',polarisation
  
  read(input_file_unit,*,END=9000,ERR=9010)fmin,fmax,fstep
  write(*,*)'fmin=',fmin
  write(*,*)'fmax=',fmax
  write(*,*)'fstep=',fstep
  write(info_file_unit,*)'fmin=',fmin
  write(info_file_unit,*)'fmax=',fmax
  write(info_file_unit,*)'fstep=',fstep
  
  call allocate_memory()
  
  do layer=1,m
  
    write(*,*)             ' '
    write(info_file_unit,*)' '
    
    write(*,*)             'reading layer number',layer
    write(info_file_unit,*)'reading layer number',layer
  
    read(input_file_unit,*,END=9000,ERR=9020),layer_number
    if (layer_number.ne.layer) then
      write(*,*)'Layers out of order, expecting layer',layer
      write(*,*)'Layer numbered ',layer_number
      stop
    end if

    read(input_file_unit,'(A1)',END=9000,ERR=9020),layer_type_in
    
    if ( (layer_type_in.eq."M").OR.(layer_type_in.eq."m") ) then
    
      layer_type(layer)=material
      write(*,*)             'Layer type=Material'
      write(info_file_unit,*)'Layer type=Material'
      
    else if ( (layer_type_in.eq."T").OR.(layer_type_in.eq."t") ) then
      layer_type(layer)=thin_layer   
      write(*,*)             'Layer type=Thin_layer'
      write(info_file_unit,*)'Layer type=Thin_layer'
      
    else
    
      write(*,*)'Unrecognised layer type'
      write(*,*)'layer type should be MATERIAL or THIN_LAYER'
      stop     
       
    end if
    
    if (layer_type(layer).eq.material) then
! read material data
      read(input_file_unit,*,END=9000,ERR=9020),layer_thickness(layer)
      write(*,*)             'Material layer thickness=',layer_thickness(layer)
      write(info_file_unit,*)'Material layer thickness=',layer_thickness(layer)
   

! read material data
      
      read(input_file_unit,'(A256)',END=9000,ERR=9020)material_filename
  
      write(*,*)'Material file name:',material_filename
      write(*,*)'Opening file:',trim(material_filename)//material_file_extn

      open(UNIT=material_file_unit,						&
           FILE=trim(material_filename)//material_file_extn,	&
           STATUS='old',							&
           ERR=9030)
      
      read(material_file_unit,*)    ! read material name
! read frequency range of validity
      read(material_file_unit,*)material_list(layer)%fmin,material_list(layer)%fmax
      read(material_file_unit,*)    ! read comment line

      call read_Sfilter(filter_in,material_file_unit) ! read permittivity filter
      material_list(layer)%eps_S=filter_in
      read(material_file_unit,*)material_list(layer)%sigma_e   
      
      read(material_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,material_file_unit) ! read permeability filter
      material_list(layer)%mu_S=filter_in                              
      read(material_file_unit,*)material_list(layer)%sigma_m   
      
      close(UNIT=material_file_unit)
      
      write(info_file_unit,*)'Material filter'
      write(info_file_unit,*)'fmin, fmax:'
      write(info_file_unit,*)material_list(layer)%fmin,material_list(layer)%fmax
      write(info_file_unit,*)'Permittivity filter:'
      call write_Sfilter(material_list(layer)%eps_S,info_file_unit)
      write(info_file_unit,*)'Electric conductivity:'
      write(info_file_unit,*)material_list(layer)%sigma_e
      write(info_file_unit,*)'Permeability filter:'
      call write_Sfilter(material_list(layer)%mu_S,info_file_unit)
      write(info_file_unit,*)'Magnetic conductivity:'
      write(info_file_unit,*)material_list(layer)%sigma_m
     
    else if (layer_type(layer).eq.thin_layer) then
    
! read thin_layer data
      layer_thickness(layer)=0d0
! read thin layer material file...    
      read(input_file_unit,'(A256)',END=9000,ERR=9020)thin_layer_filename
  
      write(*,*)'thin_layer file name:',thin_layer_filename
      write(*,*)'Opening file:',trim(thin_layer_filename)//thin_layer_file_extn

      open(UNIT=thin_layer_file_unit,						&
           FILE=trim(thin_layer_filename)//thin_layer_file_extn,	&
           STATUS='old',							&
           ERR=9050)

      read(thin_layer_file_unit,*)    ! read thin_layer name
      
! read frequency range of validity
      read(thin_layer_file_unit,*)thin_layer_list(layer)%fmin,thin_layer_list(layer)%fmax

      read(thin_layer_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,thin_layer_file_unit) ! read Z11 filter
      thin_layer_list(layer)%Z11_S=filter_in

      read(thin_layer_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,thin_layer_file_unit) ! read Z12 filter
      thin_layer_list(layer)%Z12_S=filter_in

      read(thin_layer_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,thin_layer_file_unit) ! read Z21 filter
      thin_layer_list(layer)%Z21_S=filter_in

      read(thin_layer_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,thin_layer_file_unit) ! read Z22 filter
      thin_layer_list(layer)%Z22_S=filter_in
      
      close(UNIT=thin_layer_file_unit)
      
      write(info_file_unit,*)'Thin layer filters'
      write(info_file_unit,*)'fmin, fmax:'
      write(info_file_unit,*)thin_layer_list(layer)%fmin,thin_layer_list(layer)%fmax
      write(info_file_unit,*)'Z11 filter:'
      call write_Sfilter(thin_layer_list(layer)%Z11_S,info_file_unit)
      write(info_file_unit,*)'Z12 filter:'
      call write_Sfilter(thin_layer_list(layer)%Z12_S,info_file_unit)
      write(info_file_unit,*)'Z21 filter:'
      call write_Sfilter(thin_layer_list(layer)%Z21_S,info_file_unit)
      write(info_file_unit,*)'Z22 filter:'
      call write_Sfilter(thin_layer_list(layer)%Z22_S,info_file_unit)
      
    end if  ! thin_layer
  
  end do  ! next layer
  
  write(info_file_unit,*)' '
  write(info_file_unit,*)'Finished reading input file...'
  write(info_file_unit,*)' '

  return

9000 continue
  write(*,*)'Error reading input file: End of file'
  stop

9010 continue
  write(*,*)'Error reading input file header information'
  stop

9020 continue
  write(*,*)'Error reading input file layer information, layer=',layer
  stop
  
9030 continue
  write(*,*)'Error opening material file'
  stop
  
9040 continue
  write(*,*)'Error reading material file'
  stop
  
9050 continue
  write(*,*)'Error opening thin_layer file'
  stop
  
9060 continue
  write(*,*)'Error reading thin_layer file'
  stop

  END SUBROUTINE read_Mlayer_input_file

! NAME
!     SUBROUTINE  close_files
!
! DESCRIPTION
! 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
  SUBROUTINE close_files()
  
USE Mlayer_file_module
  
! START

  write(info_file_unit,*)'Closing files'  
  
  close(UNIT=input_file_unit)       
  close(UNIT=info_file_unit)       
  close(UNIT=S11_file_unit)       
  close(UNIT=S12_file_unit)       
  close(UNIT=S21_file_unit)       
  close(UNIT=S22_file_unit)       
  close(UNIT=Z11_file_unit)       
  close(UNIT=Z12_file_unit)    
  close(UNIT=Z21_file_unit)  
  close(UNIT=Z22_file_unit)
  close(UNIT=power_file_unit)
  
  RETURN

  END SUBROUTINE close_files

