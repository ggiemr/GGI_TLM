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
! SUBROUTINE read_ngspice_output_node_list
!
! NAME
!     read_ngspice_output_node_list
!
! DESCRIPTION
!     read ngspice node output list
!
! Example packet:
!
!ngspice_node_output_list
!2   # number of ngspice output nodes
!1   # NGSPICE OUTPUT NUMBER
!1   # ngspice node number
!1   # NGSPICE OUTPUT NUMBER
!2   # ngspice node number
!
! COMMENTS
!     
!
! HISTORY
!
!     started 13/03/2019 CJS
!
!
SUBROUTINE read_ngspice_output_node_list

USE TLM_general
USE file_information
USE geometry
USE TLM_output
USE cell_parameters

IMPLICIT NONE

! local variables

integer	:: output_number
integer	:: read_number

character*256	:: input_line

! START  

  CALL write_line('CALLED: read_ngspice_output_node_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_ngspice_output_nodes
  
  CALL write_line_integer('number of output points',n_ngspice_output_nodes,0,output_to_screen_flag)
  
  if ( allocated( ngspice_output_nodes ) ) GOTO 9000
  
  allocate ( ngspice_output_nodes(1:n_ngspice_output_nodes,2) )

  do output_number=1,n_ngspice_output_nodes
  
    CALL write_line_integer('Reading ngspice output node number',output_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.output_number) goto 9010
          
      read(input_file_unit,'(A)')input_line  ! read the nodes for this ngspice output
      
      ngspice_output_nodes(output_number,1)=0
      ngspice_output_nodes(output_number,2)=0
      
      read(input_line,*,ERR=100)ngspice_output_nodes(output_number,1),  &
                                ngspice_output_nodes(output_number,2)
      GOTO 110   ! node numbers read OK
      
100   CONTINUE

! Read only a single node for this port i.e. the reference node is node zero
      read(input_line,*,ERR=9020)ngspice_output_nodes(output_number,1)

110   CONTINUE
    
        
  end do ! next output point

  CALL write_line('FINISHED: read_ngspice_output_node_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating ngspice_output_node_list:',0,.TRUE.)
     CALL write_line('ngspice_output_node_list already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading ngspice_output_node_list packet data from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
9010 CALL write_line('Error reading ngspice_output_node_list packet',0,.TRUE.)
     CALL write_line('output nodes should be numbered in order',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9020 CALL write_line('Error reading the Ngspice output node numbers',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
     
  
END SUBROUTINE read_ngspice_output_node_list
