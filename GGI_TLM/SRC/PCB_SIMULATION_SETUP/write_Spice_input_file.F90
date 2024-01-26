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
!     SUBROUTINE write_Spice_input_file
!
! DESCRIPTION
!	
!     
! COMMENTS
!     
!
!
!
! HISTORY
!
!     started 14/3/19 CJS
!     
!

SUBROUTINE write_Spice_input_file

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i

character*20 :: GGI_TLM_voltage_source_name
character*20 :: GGI_TLM_resistance_name

integer :: ngspice_reference_node
integer :: ngspice_link_node1,ngspice_link_node2
integer :: ngspice_model_internal_node

character(LEN=256) :: line

! START

  write(30,'(A)')'Ngspice template file for GGI_TLM - ngspice linked simulation'
  write(30,'(A)')'*'
  
  ngspice_reference_node=0
  
  do i=1,tot_n_ngspice_ports
  
! Each port is a link between Ngspice and GGI_TLM

    ngspice_link_node1=ngspice_port_to_node_list(i,1)
    ngspice_link_node2=ngspice_port_to_node_list(i,2)
    
    ngspice_model_internal_node=1000+i
    
    write(30,'(A20,I4,A13,I4,A4,I4)')'* GGI_TLM link port ',i,' using nodes ',ngspice_link_node1,' and',ngspice_link_node2
    
    write(GGI_TLM_voltage_source_name,'(A4,I0)')'Vtlm',i
    write(GGI_TLM_resistance_name,'(A4,I0)')'Rtlm',i
    
    write(30,'(A)')'*'
    write(30,'(A)')'* Voltage source with series resistance: equivalent circuit of TLM link'
    write(30,*)trim(GGI_TLM_voltage_source_name), ngspice_model_internal_node, ngspice_link_node2, ' DC  0.0'
    write(30,*)trim(GGI_TLM_resistance_name)    , ngspice_model_internal_node, ngspice_link_node1,      ' #Z0_TLM'
    write(30,'(A)')'* '
  
  end do ! next GGI_TLM port linked to ngspice

! Additional text from the GGI_TLM_create_PCB_simulation_model input file
  
  write(30,'(A)')'* '
  write(30,'(A)')'* Voltage source required for the voltage source controlling the breakpoint time'
  write(30,'(A)')'Vbreak time_node 0 DC 0.0'
  write(30,'(A)')'* '

! Additional text from the GGI_TLM_create_PCB_simulation_model input file

10  read(10,'(A)',END=1000,ERR=1000)line
  
    if(  index(line,'* START of Ngspice').EQ.0 ) GOTO 10    

20  read(10,'(A)',END=1000,ERR=1000)line
  
    if ( index(line,'* END of Ngspice').EQ.0 ) then
    
      write(30,'(A)')trim(line)   
      GOTO 20
      
    end if
    
1000 CONTINUE
  
  write(30,'(A)')'*'
  write(30,'(A)')'* Control for transient simulation'
!  write(30,'(A)')'.OPTIONS INTERP'
  write(30,'(A)')'.TRAN #dt_out  #tmax_ngspice 0.0 #dt_ngspice  UIC'
  write(30,'(A)')'*'
  write(30,'(A)')'.END'

RETURN  
  
END SUBROUTINE write_Spice_input_file
