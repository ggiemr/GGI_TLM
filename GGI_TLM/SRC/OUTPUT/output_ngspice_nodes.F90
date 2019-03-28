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
! SUBROUTINE initialise_ngspice_output_nodes
! SUBROUTINE write_ngspice_output_node
!
!
! NAME
!     initialise_ngspice_output_nodes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/03/2019 CJS
!
!
SUBROUTINE initialise_ngspice_output_nodes

USE TLM_general
USE TLM_output
USE ngspice_F90
USE file_information

IMPLICIT NONE

! local variables

! START
  
  if (n_ngspice_output_nodes.GT.0) then
  
    if (rank.eq.0) then
    
      OPEN(unit=ngspice_output_unit,file=trim(problem_name)//ngspice_output_extn)
  
      CALL write_time_domain_header_data(ngspice_output_unit,n_ngspice_output_nodes,n_timesteps)
      
    end if
    
  end if ! n_ngspice_output_nodes>0

  RETURN

END SUBROUTINE initialise_ngspice_output_nodes
!
! SUBROUTINE write_ngspice_output_nodes
!
! NAME
!     write_ngspice_output_nodes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/03/2019 CJS
!
!
SUBROUTINE write_ngspice_output_nodes

USE TLM_general
USE TLM_output
USE ngspice_F90
USE file_information
USE output_formats

IMPLICIT NONE

! local variables

  integer 	:: ngspice_output_node
  integer 	:: spice_node1,spice_node2
  integer 	:: opnode1,opnode2
  real*8        :: value,value1,value2

! START
  
  CALL write_line('CALLED: write_ngspice_output_nodes',0,timestepping_output_to_screen_flag)
  
  if (rank.EQ.0) then
  
    do ngspice_output_node=1,n_ngspice_output_nodes
    
       spice_node1=ngspice_output_nodes(ngspice_output_node,1)
       spice_node2=ngspice_output_nodes(ngspice_output_node,2)
       
       if ( (spice_node1.LT.0).OR.(spice_node1.GT.100) ) then
          write(*,*)'Ngspice output node1 is out of range (0-100)',spice_node1
          STOP
       end if
                
       if ( (spice_node2.LT.0).OR.(spice_node2.GT.100) ) then
          write(*,*)'Ngspice output node2 is out of range (0-100)',spice_node2
          STOP
       end if
      
       if (spice_node1.NE.0) then         
         opnode1=ngspice_node_to_V_ngspice_array_list(spice_node1)
       else
         opnode1=0
       end if
      
       if (spice_node2.NE.0) then         
         opnode2=ngspice_node_to_V_ngspice_array_list(spice_node2)
       else
         opnode2=0
       end if
                
       if ( (opnode1.LT.0).OR.(opnode1.GT.100) ) then
          write(*,*)'Ngspice output node1 is out of range (0-100)',opnode1
          STOP
       end if
                
       if ( (opnode2.LT.0).OR.(opnode2.GT.100) ) then
          write(*,*)'Ngspice output node2 is out of range (0-100)',opnode2
          STOP
       end if
       
       if (opnode1.NE.0) then
         value1=V_ngspice_array_F90(opnode1)   
       else
         value1=0d0
       end if
       
       if (opnode2.NE.0) then
         value2=V_ngspice_array_F90(opnode2)   
       else
         value2=0d0
       end if
       
       value=value1-value2
       
       if ( abs(value).lt.1D-30 )value=0d0

       write(ngspice_output_unit,time_domain_output_format)time,ngspice_output_node,value
  
    end do ! next output node
      
  end if ! (rank.eq.0)
  
  CALL write_line('FINISHED: write_ngspice_output_nodes',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE write_ngspice_output_nodes
