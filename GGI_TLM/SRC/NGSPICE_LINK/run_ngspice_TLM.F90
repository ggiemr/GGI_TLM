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
! MODULE ngspice_F90
! SUBROUTINE initialise_ngspice
! SUBROUTINE update_ngspice
!
! NAME
!     
!
! DESCRIPTION
!     
!     
! COMMENTS
!     We have hardwired the maximum number of ngpsice nodes linked to the GGI_TLM solution to be 100
!
! HISTORY
!
!     started 12/03/2019 CJS
!
MODULE ngspice_F90

USE iso_c_binding

IMPLICIT NONE

  interface
    function ngspice_wrapper_load_circuit ( ngspice_error_flag ) bind ( c )
      use iso_c_binding       
      integer ( c_int ) :: ngspice_wrapper_load_circuit
      logical (C_bool)  :: ngspice_error_flag
    end function ngspice_wrapper_load_circuit
  end interface

  interface
    function ngspice_wrapper_run ( last_time, last_V ) bind ( c )
      use iso_c_binding
      real ( c_double ) :: last_time
      real ( c_double ) :: last_V       
      integer ( c_int ) :: ngspice_wrapper_run
    end function ngspice_wrapper_run
  end interface

  interface
    function ngspice_wrapper_resume ( last_time, last_V ) bind ( c )
      use iso_c_binding
      real ( c_double ) :: last_time
      real ( c_double ) :: last_V       
      integer ( c_int ) :: ngspice_wrapper_resume
    end function ngspice_wrapper_resume
  end interface

  interface
    function ngspice_wrapper_step ( last_time, last_V, n_ngspice_nodes, ngspice_node_list, V_ngspice_array_F90 ) bind ( c )
      use iso_c_binding
      real ( c_double ) :: last_time
      real ( c_double ) :: last_V       
      integer ( c_int ) :: n_ngspice_nodes
      integer ( c_int ) :: ngspice_node_list(100)
      real ( c_double ) :: V_ngspice_array_F90(100)
      integer ( c_int ) :: ngspice_wrapper_step
    end function ngspice_wrapper_step
  end interface

  interface
    function ngspice_wrapper_run_to_breakpoint ( break_time, last_time, last_V, n_nodes, node_list, &
                                                                 V_array_F90) bind ( c )
      use iso_c_binding
      real ( c_double ) :: break_time
      real ( c_double ) :: last_time
      real ( c_double ) :: last_V       
      integer ( c_int ) :: n_nodes
      integer ( c_int ) :: node_list(100)
      real ( c_double ) :: V_array_F90(100)
      integer ( c_int ) :: ngspice_wrapper_step
   end function ngspice_wrapper_run_to_breakpoint
  end interface

  interface
    function ngspice_wrapper_set_breaktime ( breaktime ) bind ( c )
      use iso_c_binding
      real ( c_double ) :: breaktime
      integer ( c_int ) :: ngspice_wrapper_set_breaktime
    end function ngspice_wrapper_set_breaktime
  end interface
  
  INTERFACE
    INTEGER(KIND=C_INT) FUNCTION ngSpice_Command (ch_arg) BIND(C, name="ngSpice_Command")
    
      USE, INTRINSIC :: ISO_C_BINDING
      
      CHARACTER(kind=c_char), INTENT(IN) :: ch_arg(*)
      
    END FUNCTION ngSpice_Command    
  END INTERFACE

real ( c_double ) :: ngspice_break_time
real ( c_double ) :: t_ngspice_F90
real ( c_double ) :: V_ngspice_F90
integer ( c_int ) :: n_ngspice_nodes
integer ( c_int ) :: ngspice_node_list(100)
real ( c_double ) :: V_ngspice_array_F90(100)
integer :: ngspice_node_to_V_ngspice_array_list(100)
logical :: set_ngspice_node_to_V_ngspice_array_list
logical (C_bool)  :: ngspice_error_flag

real*8            :: t_eps

! Low pass filter stuff
real ( c_double ) :: V_tlm_to_ngspice_in(100,2)
real ( c_double ) :: V_tlm_to_ngspice_out(100,2)
real ( c_double ) :: V_ngspice_to_tlm_in(100,2)
real ( c_double ) :: V_ngspice_to_tlm_out(100,2)
real ( c_double ) :: LPF_k1,LPF_k2,LPF_alpha


END MODULE ngspice_F90
!
! ________________________________________________________________________________
!  
!
SUBROUTINE initialise_ngspice

! Initialise the ngspice circuit simulation which runs alongside GGI_TLM

USE constants
USE iso_c_binding
USE ngspice_F90
USE TLM_general
USE File_information

IMPLICIT NONE

! ngspice link variables
 
integer ( c_int ) :: istat

! circuit file construction stuff

character(LEN=256) :: line
character(LEN=256) :: line2
character(LEN=256) :: Z0_string
character(LEN=256) :: dt_out_string
character(LEN=256) :: dt_ngspice_string
character(LEN=256) :: tmax_string
character(LEN=256) :: command_string

logical :: Vbreak_found
  
! START

CALL write_line('CALLED: initialise_ngspice',0,output_to_screen_flag)

! set the impedance for the ngspice - GGI_TLM link. This is half the TLM link line impedance
write(Z0_string,'(ES16.6)')Z0/2d0

! set the ngspice timestep for saving data to the TLM timestep
write(dt_out_string,'(ES16.6)')dt/ngspice_timestep_factor

! set the ngspice timestep to something much less than the TLM timestep, defined by the ngspice_timestep_factor
write(dt_ngspice_string,'(ES16.6)')dt/ngspice_timestep_factor   

! Set the maximum time for the ngspice solution - add an extra bit of time to avoid problems with numerical precision
write(tmax_string,'(ES16.6)')simulation_time*1.0001d0 +dt 

! Read the template circuit file and include the parameters relating to this GGI_TLM solution
! At the same time, check that Vbreak is included and connected to time_node

open(unit=ngspice_TEMPLATE_circuit_file_unit,file='Spice_circuit_TEMPLATE.cir',status='OLD',ERR=9000)
open(unit=ngspice_circuit_file_unit,file='Spice_circuit.cir')

Vbreak_found=.FALSE.

10 CONTINUE
  read(ngspice_TEMPLATE_circuit_file_unit,'(A256)',ERR=9010,END=1000)line
    
! look for the string '#Z0_TLM' and replace with Z0
  CALL replace_in_string(line,'#Z0_TLM',Z0_string)
  
! look for the string '#dt_out' and replace with dt
  CALL replace_in_string(line,'#dt_out',dt_out_string)
  
! look for the string '#dt_ngspice' and replace with dt
  CALL replace_in_string(line,'#dt_ngspice',dt_ngspice_string)
  
! look for the string '#tmax_ngspice' and replace with tmax
  CALL replace_in_string(line,'#tmax_ngspice',tmax_string)
  
  write(ngspice_circuit_file_unit,'(A)')trim(line)
  
! convert text to lower case
  CALL convert_to_lower_case(line,256)
  if (line(1:1).NE.'*') then
    if (index(line,'vbreak').NE.0) then  
      if (index(line,'time_node').NE.0) then
        Vbreak_found=.TRUE.
      end if    
    end if    
  end if
    
  GOTO 10  ! read next line

1000 CONTINUE ! jump here when we have reached the end of the input file

close(unit=ngspice_TEMPLATE_circuit_file_unit)
close(unit=ngspice_circuit_file_unit)

if (.NOT.vbreak_found) then  

  write(*,*)'ERROR initialising Ngspice solution'
  write(*,*)'Vbreak not found in the input ciruit file'
  write(*,*)'The circuit file must include the following line to control the breakpoints:'
  write(*,*)'Vbreak time_node 0 DC 0.0'
  STOP 1
  
end if

write(*,*)' Initialise the ngspice solution'

! The c function needs to be edited to include the following:
! TLM characteristic impedance
! TLM load circuit element(s)
! Transient solution tmax and timestep (.tran) 

ngspice_error_flag=.FALSE.

istat = ngspice_wrapper_load_circuit( ngspice_error_flag )

write(*,*)'Exit status',istat
if(istat.NE.0) then
  write(*,*)'ERROR initialising Ngspice solution, istat=',istat
  STOP 1
end if
if (ngspice_error_flag) then
  write(*,*)'ERROR found in circuit file'
  STOP 1
end if

! Set the breakpoint command
write(command_string,'(A,ES16.6)')"stop when time > v(time_node)"
istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR)

! the number of nodes linking TLM and ngspice is not known until 
! the first initdata function call so set to zero for now

n_ngspice_nodes=0   

t_eps=0.5d0*dt/ngspice_timestep_factor 


LPF_alpha=0.5
LPF_k1=1.0+2.0*LPF_alpha
LPF_k2=1.0-2.0*LPF_alpha
V_tlm_to_ngspice_in(:,:)=0d0
V_tlm_to_ngspice_out(:,:)=0d0
V_ngspice_to_tlm_in(:,:)=0d0
V_ngspice_to_tlm_out(:,:)=0d0

set_ngspice_node_to_V_ngspice_array_list=.FALSE.
ngspice_node_list(:)=0
V_ngspice_array_F90(:)=0d0
ngspice_node_to_V_ngspice_array_list(:)=0

CALL write_line('FINISHED: Initialise_ngspice',0,output_to_screen_flag)

RETURN

9000 write(*,*)'ERROR: cannot open the template circuit file: Spice_circuit_TEMPLATE.cir'
STOP 1

9010 write(*,*)'ERROR reading template circuit file: Spice_circuit_TEMPLATE.cir'
STOP 1

END SUBROUTINE initialise_ngspice
!
! ________________________________________________________________________________
!  
!
SUBROUTINE update_ngspice( tlm_time )

! fortran test program to interface a 1D TLM solver with ngspice
! The 1D TLM solver simulates a transmission line between a source and a lumped circuit
! i.e. it is a terminated transmission line

USE constants
USE iso_c_binding
USE ngspice_F90
USE TLM_general
USE File_information

IMPLICIT NONE

real*8 :: tlm_time

! ngspice link variables

integer ( c_int ) :: istat

integer :: spice_node,i

character(LEN=256) :: command_string
  
! START
  
! UPDATE THE NGSPICE CIRCUIT SOLUTION
! work out the next time to halt the ngspice solution and write the breakpoint to ngspice

  ngspice_break_time=tlm_time
  
    
! write the break time to the source vbreak
  write(command_string,'(A,ES16.6)')"alter vbreak = ",ngspice_break_time
!  write(*,*)'Seeting breakpoint voltage:',trim(command_string),' nits=',ngspice_break_time/dt
  istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR)
!  write(*,*)'Returned stat:',istat
  
  istat = ngspice_wrapper_run_to_breakpoint (ngspice_break_time,t_ngspice_F90,V_ngspice_F90,              &
                                             n_ngspice_nodes,ngspice_node_list,V_ngspice_array_F90)


! ***** OLD TO BE REMOVED *****
!!  write(*,*)'CALLED update_ngspice, time=',t_ngspice_F90,' ngspice_break_time=',ngspice_break_time
!
!! single step through the ngspice solution until the time is greater than or equal to the next break time
!  do while((ngspice_break_time-t_ngspice_F90).GT.t_eps) 
!  
!    istat = ngspice_wrapper_step(t_ngspice_F90,V_ngspice_F90,n_ngspice_nodes,ngspice_node_list,V_ngspice_array_F90)
!    
!    if (istat.NE.0) then
!      write(*,*)'ERROR stepping Ngspice solution'
!      STOP 1
!    end if
!    
!  end do
! ***** END OF OLD TO BE REMOVED *****
    
  if (.NOT.set_ngspice_node_to_V_ngspice_array_list) then
! set the array which points from ngspice output array elements to V_ngspice array

    do i=1,n_ngspice_nodes  ! loop over the elements of the voltage output array
      spice_node=ngspice_node_list(i)
      
      ngspice_node_to_V_ngspice_array_list(spice_node)=i
      
    end do
    set_ngspice_node_to_V_ngspice_array_list=.TRUE.
    
  end if
  
!  write(*,*)'istat=',istat,' tout',t_ngspice_F90,' Vout',V_ngspice_F90,n_ngspice_nodes,ngspice_node_list(1),V_ngspice_array_F90(1)

! Apply a low pass filter to the voltage data going from Ngspice to TLM
    do i=1,n_ngspice_nodes  ! loop over the elements of the voltage output array

      V_ngspice_to_tlm_in(i,1)=V_ngspice_array_F90(i) 
      
! Apply first order LPF  fout(i)=(fin(i-1)+fin(i)-k2*fout(i-1))/k1

      V_ngspice_to_tlm_out(i,1)=(V_ngspice_to_tlm_in(i,2)+V_ngspice_to_tlm_in(i,1)-LPF_k2*V_ngspice_to_tlm_out(i,2))/LPF_k1
      
      V_ngspice_array_F90(i)=V_ngspice_to_tlm_out(i,1)
      
! timeshift voltage pulses
      V_ngspice_to_tlm_in(i,2)=V_ngspice_to_tlm_in(i,1)
      V_ngspice_to_tlm_out(i,2)=V_ngspice_to_tlm_out(i,1)
      
    end do


  
RETURN

END SUBROUTINE update_ngspice
