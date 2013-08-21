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
! SUBROUTINE read_volume_material_list
!
! NAME
!     read_volume_material_list
!
! DESCRIPTION
!     read volume_material list packet
!
! Example packet:
!
!
!Volume_material_list
!1   Number of volume materials (integer)
!1   VOLUME MATERIAL NUMBER
!Dispersive
!debye_1
!1       number of volumes
!1       volume list

!
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE read_volume_material_list

USE TLM_general
USE TLM_volume_materials
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_material_number
integer :: read_number
integer :: n_geometric_volumes
integer	:: i

character*256	:: input_line
character*256	:: material_name
character*256	:: material_file_name

character*256	:: material_label
  
type(SFilter)	:: filter_in

logical	:: file_exists

! function variables

character*256	:: strip_path

! START  

  CALL write_line('CALLED: Read_volume_material_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_volume_materials
  
  CALL write_line_integer('number of volume_materials',n_volume_materials,0,output_to_screen_flag)
  
  if ( allocated( volume_material_list ) ) GOTO 9000
  
  ALLOCATE( volume_material_list(1:n_volume_materials) )

  do volume_material_number=1,n_volume_materials
  
    CALL write_line_integer('Reading volume_material number',volume_material_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.volume_material_number) goto 9010
 
! read volume material type string
    read(input_file_unit,'(A)',end=9010)input_line
   
    CALL write_line('...STARTED reading volume_material_type:'//trim(input_line),0,.TRUE. )
    
! convert text to lower case
    material_name=input_line
    volume_material_list(volume_material_number)%name=material_name
    CALL convert_to_lower_case(material_name,256)
   
    if (material_name.eq.'pec') then
    
      volume_material_list(volume_material_number)%type=volume_material_type_PEC
   
    else if (material_name.eq.'pmc') then
    
      volume_material_list(volume_material_number)%type=volume_material_type_PMC
      
    else if (material_name.eq.'dispersive') then 
! dispersive thin layer, the following line of the input file should be the name of a volume material file
    
! read thin layer material file name   
      read(input_file_unit,'(A)',end=9010)input_line
      volume_material_list(volume_material_number)%name=strip_path(input_line,256)
      material_file_name=trim(input_line)//volume_material_file_extn
  
! read thin layer material file   
      CALL write_line('volume material file name=',0,output_to_screen_flag)
      CALL write_line(trim(material_file_name),0,output_to_screen_flag)

      open(UNIT=volume_material_file_unit,						&
           FILE=trim(material_file_name),	&
           STATUS='old',							&
           ERR=9020)

      volume_material_list(volume_material_number)%type=volume_material_type_DISPERSIVE
      
      read(volume_material_file_unit,'(A)')material_label    ! read material label
      
! check material label for anisotropic material
      call convert_to_lower_case(material_label,256)

! read frequency range of validity
      read(volume_material_file_unit,*)volume_material_list(volume_material_number)%fmin,	&
                                volume_material_list(volume_material_number)%fmax

      read(volume_material_file_unit,*)    ! read comment line
! read permittivity filter
      call read_Sfilter(filter_in,volume_material_file_unit) 
      volume_material_list(volume_material_number)%eps_S=filter_in
! read electric conductivity
      read(volume_material_file_unit,*)volume_material_list(volume_material_number)%sigma_e

      read(volume_material_file_unit,*)    ! read comment line
! read permeability filter
      call read_Sfilter(filter_in,volume_material_file_unit) 
      volume_material_list(volume_material_number)%mu_S=filter_in
! read magnetic conductivity
      read(volume_material_file_unit,*)volume_material_list(volume_material_number)%sigma_m
      
      close(UNIT=volume_material_file_unit)
      
      else
! not a recognised material type so an error

        GOTO 9040
	
    end if
    
! now read the volume number which have this material property
    read(input_file_unit,*,err=9030)n_geometric_volumes
    volume_material_list(volume_material_number)%n_volumes=n_geometric_volumes
    
    ALLOCATE( volume_material_list(volume_material_number)%volume_list(1:n_geometric_volumes) )
    
    read(input_file_unit,*,err=9030)( volume_material_list(volume_material_number)%volume_list(i)	&
                                     ,i=1,n_geometric_volumes )
  
    CALL write_line( '...FINISHED reading volume_material_type:'//trim(input_line),0,.TRUE. )
  
  end do ! next volume_material in volume_material_list

  CALL write_line('FINISHED: Read_volume_material_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating volume_material_list:',0,.TRUE.)
     CALL write_line('volume_material_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading volume_material_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading volume_material_list_packet_data',0,.TRUE.)
     CALL write_line('volume_materials should be numbered in order at the moment...',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading volume_material_list_packet_data',0,.TRUE.)
     CALL write_line('Problem opening volume material file:'//trim(material_file_name),0,.TRUE.)
     STOP
  
9030 CALL write_line('Error reading volume_material_list_packet_data',0,.TRUE.)
     CALL write_line('Problem reading volume material file:'//trim(material_file_name),0,.TRUE.)
     STOP
  
9040 CALL write_line('Error reading volume_material_list_packet_data',0,.TRUE.)
     CALL write_line('Unrecognised material type:'//trim(material_name),0,.TRUE.)
     STOP
    
END SUBROUTINE read_volume_material_list
