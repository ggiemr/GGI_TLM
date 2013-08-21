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
! SUBROUTINE read_cable_geometry_list
!
! NAME
!     read_cable_geometry_list
!
! DESCRIPTION
!     read cable geometry list packet. This subroutine does no processing on the cable geometry data
!     it only reads the parameters. The filling of the full cable_geometry_type structure
!     is done in SUBROUTINE create_cable_LCRG_matrices
!
! Example packet:
!
!Cable_geometry_list
!1  number of cable geometries
!1  CABLE GEOMETRY NUMBER, Cable geometry type (file name) follows:
!/home/chris/EM_MODEL_DATA/CABLE_DATA/1mm_single_wire                                       
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!
!
SUBROUTINE read_cable_geometry_list

USE TLM_general
USE Cables
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: cable_geometry_number
integer :: read_number

integer	:: n_params
integer	:: n_filters
integer :: i
  
type(SFilter)	:: filter_in

logical	:: file_exists
character(LEN=256) :: cable_geometry_filename
character(LEN=256) :: input_line

! START  

  CALL write_line('CALLED: read_cable_geometry_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_cable_geometries
  
  CALL write_line_integer('number of cable geometries',n_cable_geometries,0,output_to_screen_flag)
  
  if ( allocated( cable_geometry_list ) ) GOTO 9000
  
  allocate ( cable_geometry_list(1:n_cable_geometries) )

  do cable_geometry_number=1,n_cable_geometries
  
    CALL write_line_integer('Reading cable geometry number',cable_geometry_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.cable_geometry_number) goto 9010
 
! read cable geometry specification filename    
    read(input_file_unit,'(A)',end=9010)input_line
    
    cable_geometry_filename=trim(input_line)//cable_geometry_file_extn
   
    CALL write_line( '...STARTED reading cable geometry file:'//trim(cable_geometry_filename),0,.TRUE. )
    
    open(UNIT=cable_geometry_file_unit,FILE=cable_geometry_filename,STATUS='old',ERR=9020)

    read(cable_geometry_file_unit,'(A)',end=9030)input_line
    
! convert text to lower case
    CALL convert_to_lower_case(input_line,256)
    
    cable_geometry_list(cable_geometry_number)%cable_geometry_type_string=trim(input_line)
    CALL write_line('Reading cable type:'//trim(input_line),0,output_to_screen_flag)

    if (cable_geometry_list(cable_geometry_number)%cable_geometry_type_string.EQ.'cylindrical') then
    
      cable_geometry_list(cable_geometry_number)%cable_geometry_type=cable_geometry_type_cylindrical
      
    else if (cable_geometry_list(cable_geometry_number)%cable_geometry_type_string.EQ.'fd_cylindrical') then
    
      cable_geometry_list(cable_geometry_number)%cable_geometry_type=cable_geometry_type_FD_cylindrical
      
    else if (cable_geometry_list(cable_geometry_number)%cable_geometry_type_string.EQ.'coaxial') then
    
      cable_geometry_list(cable_geometry_number)%cable_geometry_type=cable_geometry_type_coaxial
      
    else if (cable_geometry_list(cable_geometry_number)%cable_geometry_type_string.EQ.'fd_coaxial') then
    
      cable_geometry_list(cable_geometry_number)%cable_geometry_type=cable_geometry_type_FD_coaxial
      
    else if (cable_geometry_list(cable_geometry_number)%cable_geometry_type_string.EQ.'ribbon_cable') then
      
      cable_geometry_list(cable_geometry_number)%cable_geometry_type=cable_geometry_type_ribbon
      
    else if (cable_geometry_list(cable_geometry_number)%cable_geometry_type_string.EQ.'fd_ribbon_cable') then
      
      cable_geometry_list(cable_geometry_number)%cable_geometry_type=cable_geometry_type_FD_ribbon
                
    else
    
      GOTO 9040
      
    end if

! read the number of conductors    
    read(cable_geometry_file_unit,*,err=9005)cable_geometry_list(cable_geometry_number)%n_conductors

! allocate and read parameters
    read(cable_geometry_file_unit,*,err=9005)n_params
    cable_geometry_list(cable_geometry_number)%n_parameters=n_params
    
    if (n_params.gt.0) then
    
      ALLOCATE ( cable_geometry_list(cable_geometry_number)%parameters(1:n_params) )
      
      do i=1,n_params
        read(cable_geometry_file_unit,*,err=9050)    &
                 cable_geometry_list(cable_geometry_number)%parameters(i)
      end do
      
    end if
    
! allocate and read filters
    read(cable_geometry_file_unit,*,err=9005)n_filters
    cable_geometry_list(cable_geometry_number)%n_filters=n_filters
    
    if (n_filters.gt.0) then
    
      ALLOCATE ( cable_geometry_list(cable_geometry_number)%Sfilter(1:n_filters) )
      ALLOCATE ( cable_geometry_list(cable_geometry_number)%Zfilter(1:n_filters) )
      ALLOCATE ( cable_geometry_list(cable_geometry_number)%Z_f(1:n_filters) )
      
      do i=1,n_filters
      
        read(cable_geometry_file_unit,*)    			! read comment line

        call read_Sfilter(filter_in,cable_geometry_file_unit) 	! read impedance filter
        cable_geometry_list(cable_geometry_number)%Sfilter(i)=filter_in
	
      end do
      
    end if

    close(UNIT=cable_geometry_file_unit)
      
  end do ! next cable geometry

  CALL write_line('FINISHED: read_cable_geometry_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating cable_geometry_list:',0,.TRUE.)
     CALL write_line('cable_geometry_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading cable_geometry_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading cable_geometry_list packet data',0,.TRUE.)
     CALL write_line('Cable geometries should be numbered in order',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading cable_geometry_list packet data',0,.TRUE.)
     CALL write_line('Cable geometry file des not exist. Filename:',0,.TRUE.)
     CALL write_line(trim(cable_geometry_filename),0,.TRUE.)
     STOP
  
9030 CALL write_line('Error reading cable file',0,.TRUE.)
     CALL write_line('Filename:',0,.TRUE.)
     CALL write_line(trim(cable_geometry_filename),0,.TRUE.)
     STOP
  
9040 CALL write_line('Error reading cable file',0,.TRUE.)
     CALL write_line('Filename:',0,.TRUE.)
     CALL write_line(trim(cable_geometry_filename),0,.TRUE.)
     CALL write_line('Unknown cable geometry type:'&
         //trim(cable_geometry_list(cable_geometry_number)%cable_geometry_type_string),0,.TRUE.)
     STOP
  
9050 CALL write_line('Error reading cable file',0,.TRUE.)
     CALL write_line('Filename:',0,.TRUE.)
     CALL write_line(trim(cable_geometry_filename),0,.TRUE.)
     CALL write_line('Error reading parameters',0,.TRUE.)
     STOP
  
  
END SUBROUTINE read_cable_geometry_list
