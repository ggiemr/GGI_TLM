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
! SUBROUTINE read_surface_material_list
!
! NAME
!     read_surface_material_list
!
! DESCRIPTION
!     read surface_material list packet
!
! Example packet:
!
!
!Surface_material_list
!1   Number of surface materials (integer)
!1   SURFACE NUMBER 1
!PEC
!3       number of surfaces
!1 3 4       surface list
!1 1 -1      surface orientation list
!2   SURFACE NUMBER 2
!Dispersive
!CFC
!1       number of surfaces
!2       surface list
!-1      surface orientation list
!
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE read_surface_material_list

USE TLM_general
USE TLM_surface_materials
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: surface_material_number
integer :: read_number
integer :: n_geometric_surfaces
integer	:: i

character*256	:: input_line
character*256	:: material_name
character*256	:: material_file_name

character*256	:: material_label
  
type(SFilter)	:: filter_in

! function variables

character*256	:: strip_path

logical	:: file_exists

! START  

  CALL write_line('CALLED: Read_surface_material_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_surface_materials
  
  CALL write_line_integer('number of surface_materials',n_surface_materials,0,output_to_screen_flag)
  
  if ( allocated( surface_material_list ) ) GOTO 9000
  
  ALLOCATE( surface_material_list(1:n_surface_materials) )

  do surface_material_number=1,n_surface_materials
  
    CALL write_line_integer('Reading surface_material number',surface_material_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.surface_material_number) goto 9010
 
! read surface material type string
    read(input_file_unit,'(A)',end=9010)input_line
   
    CALL write_line('...STARTED reading surface_material_type:'//trim(input_line),0,.TRUE. )
    
! convert text to lower case
    material_name=input_line
    surface_material_list(surface_material_number)%name=material_name
    CALL convert_to_lower_case(material_name,256)
   
    if (material_name.eq.'pec') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_PEC
   
    else if (material_name.eq.'pmc') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_PMC
   
    else if (material_name.eq.'free_space') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_FREE_SPACE
      
    else if (material_name.eq.'dispersive') then 
! dispersive thin layer, the following line of the input file should be the name of a surface material file
    
! read thin layer material file name   
      read(input_file_unit,'(A)',end=9010)input_line
      surface_material_list(surface_material_number)%name=strip_path(input_line,256)
      material_file_name=trim(input_line)//surface_material_file_extn
  
! read thin layer material file   
      CALL write_line('Surface material file name=',0,output_to_screen_flag)
      CALL write_line(trim(material_file_name),0,output_to_screen_flag)

      open(UNIT=surface_material_file_unit,						&
           FILE=trim(material_file_name),	&
           STATUS='old',							&
           ERR=9020)

      surface_material_list(surface_material_number)%type=surface_material_type_DISPERSIVE
      
      read(surface_material_file_unit,'(A)')material_label    ! read material label
      
! check material label for anisotropic material
      call convert_to_lower_case(material_label,256)

! read frequency range of validity
      read(surface_material_file_unit,*)surface_material_list(surface_material_number)%fmin,	&
                                surface_material_list(surface_material_number)%fmax

      read(surface_material_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z11 filter
      surface_material_list(surface_material_number)%Z11_S=filter_in

      read(surface_material_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z12 filter
      surface_material_list(surface_material_number)%Z12_S=filter_in

      read(surface_material_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z21 filter
      surface_material_list(surface_material_number)%Z21_S=filter_in

      read(surface_material_file_unit,*)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z22 filter
      surface_material_list(surface_material_number)%Z22_S=filter_in
      
      close(UNIT=surface_material_file_unit)
      
      else
! not a recognised material type so an error

        GOTO 9040
	
    end if
    
! now read the surface number which have this material property
    read(input_file_unit,*,err=9005)n_geometric_surfaces
    surface_material_list(surface_material_number)%n_surfaces=n_geometric_surfaces
    
    ALLOCATE( surface_material_list(surface_material_number)%surface_list(1:n_geometric_surfaces) )
    ALLOCATE( surface_material_list(surface_material_number)%surface_orientation_list(1:n_geometric_surfaces) )
    
    read(input_file_unit,*,err=9005)( surface_material_list(surface_material_number)%surface_list(i)	&
                                     ,i=1,n_geometric_surfaces )
    read(input_file_unit,*,err=9005)( surface_material_list(surface_material_number)%surface_orientation_list(i)	&
                                     ,i=1,n_geometric_surfaces )

! check orientation flag is +1 or -1				     
    do i=1,n_geometric_surfaces
      if (abs(surface_material_list(surface_material_number)%surface_orientation_list(i)).ne.1) goto 9050
    end do
    
    CALL write_line( '...FINISHED reading surface_material_type:'//trim(input_line),0,.TRUE. )
  
  end do ! next surface_material in surface_material_list

  CALL write_line('FINISHED: Read_surface_material_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating surface_material_list:',0,.TRUE.)
     CALL write_line('surface_material_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading surface_material_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('surface_materials should be numbered in order at the moment...',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('Problem opening surface material file:'//trim(material_file_name),0,.TRUE.)
     STOP
  
9040 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('Unrecognised material type:'//trim(material_name),0,.TRUE.)
     STOP
  
9050 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('Surface orientation should be +1 or -1',0,.TRUE.)
     STOP
    
END SUBROUTINE read_surface_material_list
