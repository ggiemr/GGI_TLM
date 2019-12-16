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
! SUBROUTINE update_TLM_to_ngspice
! SUBROUTINE update_ngspice
! SUBROUTINE update_ngspice_to_TLM
! SUBROUTINE ngspice_LPF
! SUBROUTINE finish_ngspice
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
    function ngspice_wrapper_run_to_breakpoint ( break_time, last_time, n_nodes, node_list, V_array_F90) bind ( c )
      use iso_c_binding
      real ( c_double ) :: break_time
      real ( c_double ) :: last_time
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

integer ( c_int ) :: n_ngspice_nodes
integer ( c_int ) :: ngspice_node_list(100)
real ( c_double ) :: V_ngspice_array_F90(100)

integer :: ngspice_node_to_V_ngspice_array_list(100)

logical :: set_ngspice_node_to_V_ngspice_array_list
logical (C_bool)  :: ngspice_error_flag

! lists to go from ngspice port number 
integer :: ng_material_number(100)
integer :: ng_sign(100)
integer :: ng_spice_port(100)
integer :: ng_spice_node1(100),ng_spice_node2(100)
integer :: ng_face1(100),ng_cx1(100),ng_cy1(100),ng_cz1(100)
integer :: ng_face2(100),ng_cx2(100),ng_cy2(100),ng_cz2(100)

! Checklist of ngspice port to node numbering derived from the Spice_circuit_TEMPLATE.cir file

integer :: cir_port_set(100)
integer :: cir_port_to_node1(100)
integer :: cir_port_to_node2(100)

!

real*8            :: t_eps

! Low pass filter stuff
real ( c_double ) :: V_tlm_to_ngspice_in(100,2)
real ( c_double ) :: V_tlm_to_ngspice_out(100,2)
real ( c_double ) :: V_ngspice_to_tlm_in(100,2)
real ( c_double ) :: V_ngspice_to_tlm_out(100,2)
real ( c_double ) :: LPF_k1,LPF_k2


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

integer :: n_V_ports_found
integer :: n_Z_ports_found
integer :: port_from_V(100)
integer :: node2_from_V(100)
integer :: node1_from_Z(100)
integer :: nodei_from_V(100)
integer :: nodei_from_Z(100)

logical :: Vbreak_found

integer :: i,ii,port,v1,v2,vi
logical :: Vi_Z_found  
  
! START

CALL write_line('CALLED: initialise_ngspice',0,output_to_screen_flag)

! set the ngspice timestep in relation to the TLM timestep
ngspice_dt=dt/ngspice_timestep_factor

! set the impedance for the ngspice - GGI_TLM link. This is half the TLM link line impedance
write(Z0_string,'(ES16.6)')Z0/2d0

! set the ngspice timestep for saving data to the TLM timestep
write(dt_out_string,'(ES16.6)')dt/ngspice_timestep_factor

! set the ngspice timestep to something much less than the TLM timestep, defined by the ngspice_timestep_factor
write(dt_ngspice_string,'(ES16.6)')ngspice_dt 

! Set the maximum time for the ngspice solution - add an extra bit of time to avoid problems with numerical precision
write(tmax_string,'(ES16.6)')simulation_time*1.0001d0 +dt 

! Read the template circuit file and include the parameters relating to this GGI_TLM solution
! At the same time, check that Vbreak is included and connected to time_node

n_V_ports_found=0
n_Z_ports_found=0

open(unit=ngspice_TEMPLATE_circuit_file_unit,file='Spice_circuit_TEMPLATE.cir',status='OLD',ERR=9000)
open(unit=ngspice_circuit_file_unit,file='Spice_circuit.cir')

Vbreak_found=.FALSE.

10 CONTINUE
  read(ngspice_TEMPLATE_circuit_file_unit,'(A256)',ERR=9010,END=1000)line
  
  CALL check_for_port_V(line,n_V_ports_found,port_from_V,node2_from_V,nodei_from_V)
  CALL check_for_port_Z(line,n_Z_ports_found,node1_from_Z,nodei_from_Z)
    
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

! Work out the ngspice node numbers for each port

cir_port_set(:)=0
cir_port_to_node1(:)=0
cir_port_to_node2(:)=0

write(*,*)'*********************************************'

write(*,*)'n_V_ports_found=',n_V_ports_found
write(*,*)'n_Z_ports_found=',n_Z_ports_found

if (n_V_ports_found.NE.n_Z_ports_found) then
  write(*,*)'ERROR in initialise_ngspice'
  write(*,*)'n_V_ports_found.NE.n_Z_ports_found'
  STOP 1
end if

do i=1,n_V_ports_found
  write(*,*)'Vport=',i,' port=',port_from_V(i),' node i=',nodei_from_V(i),' node 2=',node2_from_V(i)
end do
do i=1,n_Z_ports_found
  write(*,*)'Zport=',i,' node i=',nodei_from_Z(i),' node 1=',node1_from_Z(i)
end do

do i=1,n_V_ports_found

  port=port_from_V(i)
  vi=nodei_from_V(i)
  v2=node2_from_V(i)

! loop through the Z ports looking for node vi. This will allow us to determine node 2 for the port
  Vi_Z_found=.FALSE.

  do ii=1,n_Z_ports_found
    
    if (nodei_from_Z(ii).EQ.vi) then
    
      if (Vi_Z_found) then
        write(*,*)'ERROR in initialise_ngspice'
        write(*,*)'Internal port node ',vi,' found in more than one port impedance component, #z0_tlm'
        STOP 1
      end if
      
      v1=node1_from_Z(ii)
      Vi_Z_found=.TRUE.
      
    end if
      
  end do ! next Z port
  
  if (.NOT.Vi_Z_found) then
    write(*,*)'ERROR in initialise_ngspice'
    write(*,*)'Internal port node ',vi,' not found in a port impedance component, #z0_tlm'
    STOP 1  
  end if
  
  cir_port_set(port)=1
  cir_port_to_node1(port)=v1
  cir_port_to_node2(port)=v2
  
end do ! next Vport

do i=1,n_V_ports_found
  if (cir_port_set(port).EQ.1) then
    write(*,*)'port:',i,' node 1=',cir_port_to_node1(i),' node 2=',cir_port_to_node2(i)
  end if
end do

write(*,*)'*********************************************'

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

ngspice_time=0d0
last_ngspice_time=ngspice_time-ngspice_dt

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
SUBROUTINE update_TLM_to_ngspice( )

! Transfer the voltege pulses from TLM to the Ngspice model

USE constants
USE iso_c_binding
USE ngspice_F90
USE TLM_general
USE mesh
USE TLM_surface_materials
USE File_information

IMPLICIT NONE

! ngspice link variables

integer ( c_int ) :: istat

integer :: spice_port_count

integer :: sign,spice_port
integer :: node1,node2
integer :: face1,cx1,cy1,cz1
integer :: face2,cx2,cy2,cz2

real*8 :: Vspice

character(LEN=256) :: command_string
  
! START
  
! Loop over all the ports linking Ngspice to GGI_TLM

do spice_port_count=1,n_spice_ports

  sign=ng_sign(spice_port_count)
  spice_port=ng_spice_port(spice_port_count)
  node1=ng_spice_node1(spice_port_count)
  node2=ng_spice_node2(spice_port_count)
  face1=ng_face1(spice_port_count)
  cx1=ng_cx1(spice_port_count)
  cy1=ng_cy1(spice_port_count)
  cz1=ng_cz1(spice_port_count)
  face2=ng_face2(spice_port_count)
  cx2=ng_cx2(spice_port_count)
  cy2=ng_cy2(spice_port_count)
  cz2=ng_cz2(spice_port_count)
  
! Transfer the TLM incident voltage pulse(s) to the ngspice circuit

  command_string=''
  Vspice=sign*(V(face1,cx1,cy1,cz1)  + V(face2,cx2,cy2,cz2))
  if (abs(Vspice).LT.small) Vspice=0d0

! Apply a low pass filter to the voltage data going from TLM to Ngspice 
  V_tlm_to_ngspice_in(spice_port,1)=Vspice      
                 
! Apply first order LPF  fout(t)=(fin(t-1)+fin(t)-k2*fout(t-1))/k1

  CALL ngspice_LPF(V_tlm_to_ngspice_out(spice_port,1),V_tlm_to_ngspice_out(spice_port,2)  &
                  ,V_tlm_to_ngspice_in(spice_port,1),V_tlm_to_ngspice_in(spice_port,2)   )
                  
  Vspice=V_tlm_to_ngspice_out(spice_port,1)
                                       
! Build the command string including the spice node number                 
  write(command_string,'(A10,I0,A3,ES16.6)')"alter vtlm",spice_port," = ",Vspice
                 
  istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR); 
                              	      
end do ! next Spice port
  
RETURN

END SUBROUTINE update_TLM_to_ngspice
!
! ________________________________________________________________________________
!  
!
SUBROUTINE update_ngspice( tlm_time )

!  Run Ngspice for one timestep

USE constants
USE iso_c_binding
USE ngspice_F90
USE TLM_general
USE TLM_surface_materials
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
  istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR)
  
  istat = ngspice_wrapper_run_to_breakpoint (ngspice_break_time,t_ngspice_F90,n_ngspice_nodes,   &
                                             ngspice_node_list,V_ngspice_array_F90)
  
RETURN

END SUBROUTINE update_ngspice
!
! ________________________________________________________________________________
!  
!
SUBROUTINE update_ngspice_to_TLM( )

! Transfer scatterd voltages from Ngspice to TLM i.e. implement the connect process for the SPICE surface type

USE constants
USE iso_c_binding
USE ngspice_F90
USE TLM_general
USE mesh
USE TLM_surface_materials
USE File_information

IMPLICIT NONE

! ngspice link variables

integer ( c_int ) :: istat

integer :: spice_port_count

integer :: sign,spice_port
integer :: node1,node2
integer :: opnode1,opnode2
integer :: face1,cx1,cy1,cz1
integer :: face2,cx2,cy2,cz2

real*8 :: Vspice,Vspice1,Vspice2

integer :: spice_node

character(LEN=256) :: command_string
integer :: i
  
! START

!write(*,*)'CALLED: update_ngspice_to_TLM'

! Wait until the Ngspice solution has reached the breakpoint then transfer node voltages to V_ngspice_array_F90
  
  istat = ngspice_wrapper_run_to_breakpoint (ngspice_break_time,t_ngspice_F90,n_ngspice_nodes,   &
                                             ngspice_node_list,V_ngspice_array_F90) 

! Note that Ngspice has to be run at least once before we can carry out the following process
! Therefore it needs to be put here after ngspice_wrapper_run_to_breakpoint    

  if (.NOT.set_ngspice_node_to_V_ngspice_array_list) then
  
! set the array which points from ngspice output array elements to V_ngspice array

    do i=1,n_ngspice_nodes  ! loop over the elements of the voltage output array
      spice_node=ngspice_node_list(i)
      
      ngspice_node_to_V_ngspice_array_list(spice_node)=i
      
    end do
    set_ngspice_node_to_V_ngspice_array_list=.TRUE.
    
  end if

! Apply a low pass filter to the voltage data going from Ngspice to TLM

    do i=1,n_ngspice_nodes  ! loop over the elements of the voltage output array

      V_ngspice_to_tlm_in(i,1)=V_ngspice_array_F90(i) 
      
! Apply first order LPF  fout(i)=(fin(i-1)+fin(i)-k2*fout(i-1))/k1

      CALL ngspice_LPF(V_ngspice_to_tlm_out(i,1),V_ngspice_to_tlm_out(i,2),V_ngspice_to_tlm_in(i,1),V_ngspice_to_tlm_in(i,2))
      
    end do  
  
RETURN

END SUBROUTINE update_ngspice_to_TLM
!
! ________________________________________________________________________________
!  
!
SUBROUTINE ngspice_LPF( fout,fout_m1,fin,fin_m1  )

USE ngspice_F90

IMPLICIT NONE

real ( c_double ) :: fout,fout_m1,fin,fin_m1

! Apply a first order low pass filter to time sampled data
!  fout(i)=(fin(i-1)+fin(i)-k2*fout(i-1))/k1

! START

fout=(fin_m1+fin-LPF_k2*fout_m1)/LPF_k1

! timeshift the input and output time samples

fin_m1=fin
fout_m1=fout

RETURN

END SUBROUTINE ngspice_LPF
!
! ________________________________________________________________________________
!  
!
SUBROUTINE finish_ngspice

! Finalise the ngspice circuit simulation process which runs alongside GGI_TLM

USE constants
USE iso_c_binding
USE ngspice_F90
USE TLM_general
USE File_information

IMPLICIT NONE

integer ( c_int ) :: istat

character(LEN=256) :: command_string

! START

  write(command_string,'(A,ES16.6)')"bg_halt"
  istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR)

RETURN

END SUBROUTINE finish_ngspice
!
! ________________________________________________________________________________
!  
!
SUBROUTINE check_for_port_V(line,nports,port,node1,nodei)

IMPLICIT NONE

character(LEN=256) :: line
integer :: nports,port(100),node1(100),nodei(100)

! local variables

character(LEN=256) :: lcstring
character(LEN=256) :: substring
integer :: pos
integer :: length

! START

!write(*,*)'CALLED: check_for_port_V, nports=',nports

lcstring=line
CALL convert_to_lower_case(lcstring,256)

pos=INDEX(lcstring, 'vtlm')

if (pos.EQ.0) RETURN

nports=nports+1

if (nports.GT.100) then
  write(*,*)'ERROR in initialise_ngspice'
  write(*,*)'The number of ports linking GGI_TLM and Ngspice should be less than 100'
  write(*,*)'nports=',nports
  STOP 1
end if

!write(*,*)'**** Found GGI_TLM V port in string ****'
!write(*,*)'nports=',nports
!write(*,*)trim(lcstring)
!write(*,*)'position=',pos

length=LEN(lcstring)
substring=trim(lcstring(pos+4:length))

!write(*,*)'Checking substring for port, nodei, node1'
!write(*,*)trim(substring)

read(substring,*)port(nports), nodei(nports), node1(nports)
!write(*,*)'port=',port(nports)
!write(*,*)'nodei=',nodei(nports)
!write(*,*)'node1=',node1(nports)

END SUBROUTINE check_for_port_V
!
! ________________________________________________________________________________
!  
!
SUBROUTINE check_for_port_Z(line,nports,node2,nodei)

IMPLICIT NONE

character(LEN=256) :: line
integer :: nports,node2(100),nodei(100)

! local variables

character(LEN=256) :: lcstring
character(LEN=256) :: substring
character(LEN=256) :: subsubstring
integer :: pos
integer :: length

! START

!write(*,*)'CALLED: check_for_port_Z, nports=',nports

lcstring=line
CALL convert_to_lower_case(lcstring,256)

pos=INDEX(lcstring, '#z0_tlm')

if (pos.EQ.0) RETURN

nports=nports+1

if (nports.GT.100) then
  write(*,*)'ERROR in initialise_ngspice'
  write(*,*)'The number of ports linking GGI_TLM and Ngspice should be less than 100'
  write(*,*)'nports=',nports
  STOP 1
end if

!write(*,*)'**** Found GGI_TLM Z port in string ****'
!write(*,*)'nports=',nports
!write(*,*)trim(lcstring)
!write(*,*)'position=',pos

! find the end of the resistive element string, starting with 'r'

length=LEN(lcstring)
pos=INDEX(lcstring,'r')

if (pos.EQ.0) then
  write(*,*)'ERROR in initialise_ngspice'
  write(*,*)'Could not find port resistance specified. String:'
  write(*,*)trim(lcstring)
  STOP 1
end if

substring=lcstring(pos:length)
pos=INDEX(substring,' ')

if (pos.EQ.0) then
  write(*,*)'ERROR in initialise_ngspice'
  write(*,*)'Could not identify port resistance node numbers. String:'
  write(*,*)trim(substring)
  STOP 1
end if

length=LEN(substring)
subsubstring=substring(pos+1:length)

!write(*,*)'Checking substring for nodei, node1'
!write(*,*)trim(subsubstring)

read(subsubstring,*) nodei(nports), node2(nports)
!write(*,*)'nodei=',nodei(nports)
!write(*,*)'node2=',node2(nports)

END SUBROUTINE check_for_port_Z
