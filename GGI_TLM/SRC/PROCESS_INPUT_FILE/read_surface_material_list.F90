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
!5   Number of surface materials (integer)
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
!3       SURFACE MATERIAL NUMBER 
!DIODE
!1e-12  0.026   Diode Is and nVt values
!+x              Diode forward direction, infinite imedance in the orthogonal direction
!1       number of surfaces
!7       surface list
!1       surface orientation list
!4       SURFACE MATERIAL NUMBER 
!SPICE
!1            # ngspice circuit node number in the Spice_circuit_TEMPLATE.cir file for the port (note, assume port voltage is referenced to node zero)
!-x           # port direction, infinite impedance in the orthogonal direction
!1       number of surfaces
!7       surface list
!1       surface orientation list
!5       SURFACE MATERIAL NUMBER 
!SWITCH
!0.5e-7 0.75e-7 1e-7                                   ! time_to switch_on  time_to_switch_off period
!1       number of surfaces
!7       surface list
!1       surface orientation list
!
! COMMENTS
!     
!
! HISTORY
!
!    started 14/08/2012 CJS
!    2/12/2013 		CJS: Implement anisotropic impedance boundary conditions
!    3/09/2014		CJS: Implement simple diode impedance boundary conditions
!    11/03/2019		CJS: Implement SPICE circuit model link
!    20/03/2024		CJS: Implement SWITCH model
!    6/09/2024		CJS: Implement DSWITCH model
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

  read(input_file_unit,*,err=9005,end=9005)n_surface_materials
  
  CALL write_line_integer('number of surface_materials',n_surface_materials,0,output_to_screen_flag)
  
  if ( allocated( surface_material_list ) ) GOTO 9000
  
  ALLOCATE( surface_material_list(1:n_surface_materials) )
  
  n_spice_ports=0

  do surface_material_number=1,n_surface_materials
  
    CALL write_line_integer('Reading surface_material number',surface_material_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005,end=9005)read_number
    if (read_number.ne.surface_material_number) goto 9010
 
! read surface material type string
    read(input_file_unit,'(A)',err=9010,end=9010)input_line
   
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
! DISPERSIVE THIN LAYER, the following line of the input file should be the name of a surface material file
    
! read thin layer material file name   
      read(input_file_unit,'(A)',err=9010,end=9010)input_line
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
      
      read(surface_material_file_unit,'(A)',err=9030,end=9030)material_label    ! read material label
      
! check material label for anisotropic material
      call convert_to_lower_case(material_label,256)

! read frequency range of validity
      read(surface_material_file_unit,*,err=9030,end=9030)surface_material_list(surface_material_number)%fmin,	&
                                surface_material_list(surface_material_number)%fmax
! READ X POLARISATION FILTER
      write(*,*)'Read x polarisation filter:'

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z11 filter
      surface_material_list(surface_material_number)%Z11_S(1)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z12 filter
      surface_material_list(surface_material_number)%Z12_S(1)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z21 filter
      surface_material_list(surface_material_number)%Z21_S(1)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z22 filter
      surface_material_list(surface_material_number)%Z22_S(1)=filter_in

! COPY X POLARISATION FILTER TO Y AND Z POLARISATIONS

      surface_material_list(surface_material_number)%Z11_S(2)=surface_material_list(surface_material_number)%Z11_S(1)
      surface_material_list(surface_material_number)%Z12_S(2)=surface_material_list(surface_material_number)%Z12_S(1)
      surface_material_list(surface_material_number)%Z21_S(2)=surface_material_list(surface_material_number)%Z21_S(1)
      surface_material_list(surface_material_number)%Z22_S(2)=surface_material_list(surface_material_number)%Z22_S(1)

      surface_material_list(surface_material_number)%Z11_S(3)=surface_material_list(surface_material_number)%Z11_S(1)
      surface_material_list(surface_material_number)%Z12_S(3)=surface_material_list(surface_material_number)%Z12_S(1)
      surface_material_list(surface_material_number)%Z21_S(3)=surface_material_list(surface_material_number)%Z21_S(1)
      surface_material_list(surface_material_number)%Z22_S(3)=surface_material_list(surface_material_number)%Z22_S(1)
     
      close(UNIT=surface_material_file_unit)
      
      
    else if (material_name.eq.'anisotropic_dispersive') then 
! ANISOTROPIC_DISPERSIVE THIN LAYER, the following line of the input file should be the name of a surface material file
    
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

      surface_material_list(surface_material_number)%type=surface_material_type_ANISOTROPIC_DISPERSIVE
      
      read(surface_material_file_unit,'(A)',err=9030,end=9030)material_label    ! read material label
      
! check material label for anisotropic material
      call convert_to_lower_case(material_label,256)

! read frequency range of validity
      read(surface_material_file_unit,*)surface_material_list(surface_material_number)%fmin,	&
                                surface_material_list(surface_material_number)%fmax

! READ X POLARISATION FILTER
      write(*,*)'Read x polarisation filter:'
      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z11 filter
      surface_material_list(surface_material_number)%Z11_S(1)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z12 filter
      surface_material_list(surface_material_number)%Z12_S(1)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z21 filter
      surface_material_list(surface_material_number)%Z21_S(1)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z22 filter
      surface_material_list(surface_material_number)%Z22_S(1)=filter_in

! READ Y POLARISATION FILTER
      write(*,*)'Read y polarisation filter:'
      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z11 filter
      surface_material_list(surface_material_number)%Z11_S(2)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z12 filter
      surface_material_list(surface_material_number)%Z12_S(2)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z21 filter
      surface_material_list(surface_material_number)%Z21_S(2)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z22 filter
      surface_material_list(surface_material_number)%Z22_S(2)=filter_in

! READ Z POLARISATION FILTER
      write(*,*)'Read z polarisation filter:'
      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z11 filter
      surface_material_list(surface_material_number)%Z11_S(3)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z12 filter
      surface_material_list(surface_material_number)%Z12_S(3)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z21 filter
      surface_material_list(surface_material_number)%Z21_S(3)=filter_in

      read(surface_material_file_unit,*,err=9030,end=9030)    ! read comment line
      call read_Sfilter(filter_in,surface_material_file_unit) ! read Z22 filter
      surface_material_list(surface_material_number)%Z22_S(3)=filter_in
      
      close(UNIT=surface_material_file_unit)
   
    else if (material_name.eq.'diode') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_DIODE

! read diode parameters      
     read(input_file_unit,*,err=9005,end=9005)surface_material_list(surface_material_number)%Diode_Is,	&
                                     surface_material_list(surface_material_number)%Diode_nVt,	&
                                     surface_material_list(surface_material_number)%Diode_Rs,	&
				     surface_material_list(surface_material_number)%Diode_Cj
				     
     read(input_file_unit,'(A2)',err=9005)surface_material_list(surface_material_number)%Diode_direction
     
! check direction is OK
     if ( (surface_material_list(surface_material_number)%Diode_direction.NE.'-x').AND.	&
          (surface_material_list(surface_material_number)%Diode_direction.NE.'+x').AND. &
          (surface_material_list(surface_material_number)%Diode_direction.NE.'-y').AND.	&
          (surface_material_list(surface_material_number)%Diode_direction.NE.'+y').AND. &
          (surface_material_list(surface_material_number)%Diode_direction.NE.'-z').AND.	&
          (surface_material_list(surface_material_number)%Diode_direction.NE.'+z') ) then
	  
	write(*,*)'Diode direction shoule be -x, +x, -y, +y, -z or +z'
	STOP
	  	  
      end if
      
      surface_material_list(surface_material_number)%diode_sign=-1
      if (surface_material_list(surface_material_number)%Diode_direction(1:1).EQ.'-') then
        surface_material_list(surface_material_number)%diode_sign=+1
      end if  
   
    else if (material_name.eq.'spice') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_SPICE

! read spice circuit model parameters      

! read the Ngspice node(s) for this surface 
      n_spice_ports=n_spice_ports+1
      read(input_file_unit,*,err=9005,end=9005)surface_material_list(surface_material_number)%Spice_circuit_file_port
      
      if (surface_material_list(surface_material_number)%Spice_circuit_file_port.NE.n_spice_ports) then
        write(*,*)'ERROR in read_surface_material_list: Spice ports should be numbered in order at the moment'
        STOP 1
      end if
      
      read(input_file_unit,'(A)',err=9005,end=9005)input_line  ! read the nodes for this port
      
      surface_material_list(surface_material_number)%Spice_circuit_file_nodes(1)=0
      surface_material_list(surface_material_number)%Spice_circuit_file_nodes(2)=0
      
      read(input_line,*,ERR=100)surface_material_list(surface_material_number)%Spice_circuit_file_nodes(1),  &
                                surface_material_list(surface_material_number)%Spice_circuit_file_nodes(2)

! check that the node number is less than 100                                
      if (surface_material_list(surface_material_number)%Spice_circuit_file_nodes(1).GT.99) then
        write(*,*)'ERROR reading surface material number',surface_material_number
        write(*,*)'Spice node number should be less than 100'
        write(*,*)'Spice node number=',surface_material_list(surface_material_number)%Spice_circuit_file_nodes(1)
        STOP 1
      end if
      
      if (surface_material_list(surface_material_number)%Spice_circuit_file_nodes(2).GT.99) then
        write(*,*)'ERROR reading surface material number',surface_material_number
        write(*,*)'Spice node number should be less than 100'
        write(*,*)'Spice node number=',surface_material_list(surface_material_number)%Spice_circuit_file_nodes(2)
        STOP 1
      end if
      
      GOTO 110   ! node numbers read OK
      
100   CONTINUE

! Read only a single node for this port i.e. the reference node is node zero
      read(input_line,*,ERR=9060)surface_material_list(surface_material_number)%Spice_circuit_file_nodes(1)

! check that the node number is less than 100                                
      if (surface_material_list(surface_material_number)%Spice_circuit_file_nodes(1).GT.99) then
        write(*,*)'ERROR reading surface material number',surface_material_number
        write(*,*)'Spice node number should be less than 100'
        write(*,*)'Spice node number=',surface_material_list(surface_material_number)%Spice_circuit_file_nodes(1)
        STOP 1
      end if

110   CONTINUE
    		     
      read(input_file_unit,'(A2)',err=9005)surface_material_list(surface_material_number)%Spice_port_direction
     
! check direction is OK
      if ( (surface_material_list(surface_material_number)%Spice_port_direction.NE.'-x').AND.	&
           (surface_material_list(surface_material_number)%Spice_port_direction.NE.'+x').AND. &
           (surface_material_list(surface_material_number)%Spice_port_direction.NE.'-y').AND.	&
           (surface_material_list(surface_material_number)%Spice_port_direction.NE.'+y').AND. &
           (surface_material_list(surface_material_number)%Spice_port_direction.NE.'-z').AND.	&
           (surface_material_list(surface_material_number)%Spice_port_direction.NE.'+z') ) then
	  
	 write(*,*)'Spice_port direction shoule be -x, +x, -y, +y, -z or +z'
	 STOP
	  
       end if
      
      surface_material_list(surface_material_number)%Spice_port_sign=-1
      if (surface_material_list(surface_material_number)%Spice_port_direction(1:1).EQ.'-') then
        surface_material_list(surface_material_number)%Spice_port_sign=1
      end if  
      
      run_ngspice=.TRUE.
   
    else if (material_name.eq.'switch') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_SWITCH
      
! Read switch_t_on switch_t_off switch_period
      read(input_file_unit,*,err=9005,end=9005)surface_material_list(surface_material_number)%switch_t_on, &
                                     surface_material_list(surface_material_number)%switch_t_off,	   &
                                     surface_material_list(surface_material_number)%switch_period

      surface_material_list(surface_material_number)%switch_delay =0d0
      surface_material_list(surface_material_number)%switch_R_on =0d0
      surface_material_list(surface_material_number)%switch_R_off=1d12
      
    else if (material_name.eq.'dswitch') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_SWITCH
      
! Read switch_t_on switch_t_off switch_period and switch_delay
      read(input_file_unit,*,err=9005,end=9005)surface_material_list(surface_material_number)%switch_t_on, &
                                     surface_material_list(surface_material_number)%switch_t_off,	   &
                                     surface_material_list(surface_material_number)%switch_period,         &
                                     surface_material_list(surface_material_number)%switch_delay

      surface_material_list(surface_material_number)%switch_R_on =0d0
      surface_material_list(surface_material_number)%switch_R_off=1d12
    
    else if (material_name.eq.'rswitch') then
    
      surface_material_list(surface_material_number)%type=surface_material_type_SWITCH
      
! Read switch_t_on switch_t_off switch_period switch_R_on switch_R_off
      read(input_file_unit,*,err=9005,end=9005)surface_material_list(surface_material_number)%switch_t_on, &
                                     surface_material_list(surface_material_number)%switch_t_off,	   &
                                     surface_material_list(surface_material_number)%switch_period,	   &
                                     surface_material_list(surface_material_number)%switch_R_on,	   &
                                     surface_material_list(surface_material_number)%switch_R_off

      surface_material_list(surface_material_number)%switch_delay =0d0
      
    else
! not a recognised material type so an error

      GOTO 9040
	
    end if
    
! now read the surface number which have this material property
    read(input_file_unit,*,err=9005)n_geometric_surfaces
    surface_material_list(surface_material_number)%n_surfaces=n_geometric_surfaces
    
    ALLOCATE( surface_material_list(surface_material_number)%surface_list(1:n_geometric_surfaces) )
    ALLOCATE( surface_material_list(surface_material_number)%surface_orientation_list(1:n_geometric_surfaces) )
    
    if ( (surface_material_list(surface_material_number)%type.EQ.surface_material_type_SPICE).AND.   &
         (n_geometric_surfaces.NE.1) ) then
      
      write(*,*)'ERROR reading surface material number',surface_material_number
      write(*,*)'Spice surfaces shoould only consist of a single geometric surface'
      STOP 1
    
    end if 
    
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
     STOP 1
  
9005 CALL write_line('Error reading surface_material_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1

9010 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('surface_materials should be numbered in order at the moment...',0,.TRUE.)
     STOP 1
  
9020 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('Problem opening surface material file:'//trim(material_file_name),0,.TRUE.)
     STOP 1
    
9030 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('Error reading surface_material file',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1

9040 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('Unrecognised material type:'//trim(material_name),0,.TRUE.)
     STOP 1
  
9050 CALL write_line('Error reading surface_material_list_packet_data',0,.TRUE.)
     CALL write_line('Surface orientation should be +1 or -1',0,.TRUE.)
     STOP 1

9060 CALL write_line('Error reading the Ngspice node numbers for the surface',0,.TRUE.)
     STOP 1
    
END SUBROUTINE read_surface_material_list
