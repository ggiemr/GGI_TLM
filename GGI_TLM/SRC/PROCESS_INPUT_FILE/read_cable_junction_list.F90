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
! SUBROUTINE read_cable_junction_list
!
! NAME
!     read_cable_junction_list
!
! DESCRIPTION
!     read cable list packet
!
! Example packet:
!
!2    number of junctions
!1    JUNCTION NUMBER
!2    junction point number
!1   number of internal connection nodes (n_int)
!2   number of cables
!1 2   cable list
!2 1   corresponding cable end number list
!1   	Cable 1 number of external conductors (n_ext)
!1   	Cable 1 P matrix: matrix (n_int rows* n_ext columns)
!1   	Cable 1 voltage source function list
!25.0 	Cable 1 impedance list
!1   	Cable 2 number of external conductors (n_ext)
!1   	Cable 2 P matrix: matrix (n_int rows* n_ext columns)
!1   	Cable 2 voltage source function list
!25.0 	Cable 2 impedance list
!0	number of internal impedances
!2    JUNCTION NUMBER
!3    junction point number
!2   number of internal connection nodes (n_int)
!2   number of cables
!1 3   cable list
!2 1   corresponding cable end number list
!1   	Cable 1 number of external conductors (n_ext)
!1   	Cable 1 P matrix: matrix (n_int rows* n_ext columns)
!0 
!1   	Cable 1 voltage source function list
!5.0 	Cable 1 impedance list
!1   	Cable 2 number of external conductors (n_ext)
!0  	Cable 2 P matrix: matrix (n_int rows* n_ext columns)
!1 
!1   	Cable 2 voltage source function list
!0.0 	Cable 2 impedance list
!1   number of internal impedances
!1   Internal impedance number
!1 2 ! between internal conection nodes
!# Impedance filter function
!1.0e0 		# w normalisation constant
!1   		# a order
!0.0987 6.59e-10	# a coefficients
!0   		# b order
!1.0 		# b coefficients
!FACE
!-1  boundary condition (-1 indicates V=0 i.e. connection to a surface)
!
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS - may change packet format to include more data such as
!                              frequency dependent loads, current sources, etc
!	19/11/2012 CJS Include termination to surfaces
!       13/12/2012 CJS Start to include frequency dependent impedances at junctions
!
SUBROUTINE read_cable_junction_list

USE TLM_general
USE Cables
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: cable_junction_number
integer	:: cable_number
integer :: read_number

integer :: number_of_cables
integer :: number_of_lines
integer :: n_int,n_ext
integer :: n_Z
integer :: impedance_number
integer :: row

integer :: i

integer :: Nbc
  
type(SFilter)	:: filter_in

character*256	:: input_line
character	:: ch

! START  

  CALL write_line('CALLED: read_cable_junction_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_cable_junctions
  
  CALL write_line_integer('number of cable junctions',n_cable_junctions,0,output_to_screen_flag)
  
  if ( allocated( cable_junction_list ) ) GOTO 9005
  
  ALLOCATE ( cable_junction_list(1:n_cable_junctions) )

  do cable_junction_number=1,n_cable_junctions
  
    CALL write_line_integer('Reading cable junction number',cable_junction_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.cable_junction_number) goto 9010
 
! read the junction point number 
    read(input_file_unit,*,err=9005)cable_junction_list(cable_junction_number)%point_number
 
! read the number of internal connection nodes
    read(input_file_unit,*,err=9005)n_int
    cable_junction_list(cable_junction_number)%n_internal_connection_nodes=n_int

! read the number of cables connecting to the junction
    read(input_file_unit,*,err=9005)number_of_cables
    cable_junction_list(cable_junction_number)%number_of_cables=number_of_cables

    if (number_of_cables.gt.0) then
! allocate and read the cable data
    
      ALLOCATE ( cable_junction_list(cable_junction_number)%cable_list(1:number_of_cables) )
      ALLOCATE ( cable_junction_list(cable_junction_number)%cable_end_list(1:number_of_cables) )
      ALLOCATE ( cable_junction_list(cable_junction_number)%n_external_conductors(1:number_of_cables) )
      ALLOCATE ( cable_junction_list(cable_junction_number)%Pmatrix(1:number_of_cables) )
      ALLOCATE ( cable_junction_list(cable_junction_number)%excitation_function(1:number_of_cables) )
      ALLOCATE ( cable_junction_list(cable_junction_number)%resistance(1:number_of_cables) )

! read the cable list
      read(input_file_unit,*,err=9020)(cable_junction_list(cable_junction_number)%cable_list(i),i=1,number_of_cables)

! read the cable end list
      read(input_file_unit,*,err=9030)(cable_junction_list(cable_junction_number)%cable_end_list(i),i=1,number_of_cables)
      
      do cable_number=1,number_of_cables
      
        read(input_file_unit,*,err=9030)n_ext
	cable_junction_list(cable_junction_number)%n_external_conductors(cable_number)=n_ext
                                                   
! allocate and read P matrix
        ALLOCATE ( cable_junction_list(cable_junction_number)%Pmatrix(cable_number)%P(n_int,n_ext) )
        do row=1,n_int
	  read(input_file_unit,*,err=9040)	&
	       (cable_junction_list(cable_junction_number)%Pmatrix(cable_number)%P(row,i),i=1,n_ext)
	end do 
	
        ALLOCATE ( cable_junction_list(cable_junction_number)%excitation_function(cable_number)%value(1:n_ext) )
	read(input_file_unit,*,err=9050)	&
	          (cable_junction_list(cable_junction_number)%excitation_function(cable_number)%value(i),i=1,n_ext)
	
        ALLOCATE ( cable_junction_list(cable_junction_number)%resistance(cable_number)%value(1:n_ext) )	       
	read(input_file_unit,*,err=9060)	&
	          (cable_junction_list(cable_junction_number)%resistance(cable_number)%value(i),i=1,n_ext)
      
      end do ! next cable connected to the junction
      
    end if   ! number_of_cables.GT.0

! Read the number of internal impedances and the associated data

    read(input_file_unit,*,err=9005)n_Z
    cable_junction_list(cable_junction_number)%number_of_internal_impedances=n_Z
    if (nz.ne.0) then
    
      ALLOCATE( cable_junction_list(cable_junction_number)%node_1(1:n_Z) )
      ALLOCATE( cable_junction_list(cable_junction_number)%node_2(1:n_Z) )
      ALLOCATE( cable_junction_list(cable_junction_number)%Sfilter(1:n_Z) )
    
      do i=1,n_Z
      
        read(input_file_unit,*)impedance_number
	
	if (impedance_number.NE.i) GOTO 9080
	
        read(input_file_unit,*)cable_junction_list(cable_junction_number)%node_1(i),	&
	                       cable_junction_list(cable_junction_number)%node_2(i)
			       
        read(input_file_unit,*)			! comment line before filter data   
 
        call read_Sfilter(filter_in,input_file_unit) 	! read impedance filter
        cable_junction_list(cable_junction_number)%Sfilter(i)=filter_in
    
      end do ! next internal impedance
      
    end if ! number_of_internal_impedances.NE.0
    
! read the next line to see whether this is a face junction or a cell centre junction
! assume we have a cell centre junction initially

    cable_junction_list(cable_junction_number)%junction_type=junction_type_cell
    
    read(input_file_unit,'(A1)',err=1000)ch
    
    if ( (ch.eq.'f').OR.(ch.eq.'F') ) then
    
      cable_junction_list(cable_junction_number)%junction_type=junction_type_face
      
! allocate and read boundary condition data
      ALLOCATE ( cable_junction_list(cable_junction_number)%BC(1:n_int) )	       
      read(input_file_unit,*,err=9070)	&
	          (cable_junction_list(cable_junction_number)%BC(i),i=1,n_int)
		  
! Check that only one internal connection node is set to -1
      Nbc=0
      do i=1,n_int
	if (cable_junction_list(cable_junction_number)%BC(i).EQ.-1) then
	  Nbc=Nbc+1
	else if (cable_junction_list(cable_junction_number)%BC(i).NE.0) then
	  GOTO 9090
	end if
      end do
      
      if (Nbc.gt.1) GOTO 9100
        
    else
! this is not a face junction so go back as if the extra line had not been read
  
      backspace(unit=input_file_unit)
  
    end if
    
1000 CONTINUE
      
  end do ! next cable junction

  CALL write_line('FINISHED: read_cable_junction_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating cable_junction_list:',0,.TRUE.)
     CALL write_line('cable_junction_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading cable_junction_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line('Cable junctions should be numbered in order',0,.TRUE.)
     STOP

9020 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line_integer('Error in cable line list, cable junction:',cable_junction_number,0,.TRUE.)
     STOP

9030 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line_integer('Error in cable line end list, cable junction:',cable_junction_number,0,.TRUE.)
     STOP

9040 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line_integer('Error reading P matrix, cable junction:',cable_junction_number,0,.TRUE.)
     CALL write_line_integer('Connecting cable number:',cable_number,0,.TRUE.)
     STOP

9050 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line_integer('Error reading excitation function vector, cable junction:',cable_junction_number,0,.TRUE.)
     CALL write_line_integer('Connecting cable number:',cable_number,0,.TRUE.)
     STOP

9060 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line_integer('Error reading resistance vector, cable junction:',cable_junction_number,0,.TRUE.)
     CALL write_line_integer('Connecting cable number:',cable_number,0,.TRUE.)
     STOP

9070 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line_integer('Error reading boundary condition vector, cable junction:',cable_junction_number,0,.TRUE.)
     STOP

9080 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line('Internal impedances should be numbered in order',0,.TRUE.)
     STOP

9090 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line('Boundary conditions on internal connection nodes should be 0 or -1',0,.TRUE.)
     STOP

9100 CALL write_line('Error reading cable_junction_list packet data',0,.TRUE.)
     CALL write_line('At most only one boundary condition on internal connection nodes should be -1',0,.TRUE.)
     STOP

  
  
END SUBROUTINE read_cable_junction_list
