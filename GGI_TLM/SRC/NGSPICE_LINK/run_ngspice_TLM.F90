MODULE ngspice_F90

USE iso_c_binding

IMPLICIT NONE

  interface
    function ngspice_wrapper_load_circuit (  ) bind ( c )
      use iso_c_binding       
      integer ( c_int ) :: ngspice_wrapper_load_circuit
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
    function ngspice_wrapper_step ( last_time, last_V, n_nodes, node_list, V_array_F90 ) bind ( c )
      use iso_c_binding
      real ( c_double ) :: last_time
      real ( c_double ) :: last_V       
      integer ( c_int ) :: n_nodes
      integer ( c_int ) :: node_list(100)
      real ( c_double ) :: V_array_F90(100)
      integer ( c_int ) :: ngspice_wrapper_step
    end function ngspice_wrapper_step
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

END MODULE ngspice_F90
!
! ________________________________________________________________________________
!  
!
SUBROUTINE run_ngspice_TLM

! fortran test program to interface a 1D TLM solver with ngspice
! The 1D TLM solver simulates a transmission line between a source and a lumped circuit
! i.e. it is a terminated transmission line

USE iso_c_binding
USE ngspice_F90

IMPLICIT NONE

! ngspice link variables

integer ( c_int ) :: istat

real ( c_double ) :: break_time
real ( c_double ) :: t_F90
real ( c_double ) :: V_F90
integer ( c_int ) :: n_nodes
integer ( c_int ) :: node_list(100)
real ( c_double ) :: V_array_F90(100)
INTEGER(kind=C_INT) :: ret
character*80 :: command_string

! timestepping
real ( c_double ) :: dt,tmax,eps,time
integer :: nt,i,iz,nodeloop,opnode

! circuit file construction stuff

character(LEN=256) :: line
character(LEN=256) :: Z0_string
character(LEN=256) :: dt_string
character(LEN=256) :: tmax_string

! local TLM stuff

real( c_double )  :: dl    ! TLM cell dimension
integer :: nz    ! Number of TLM cells
real( c_double )  :: Z0    ! TLM transmission line impedance
real( c_double )  :: c0    ! TLM transmission line velocity

real( c_double )  :: Zs    ! source impedance
real( c_double )  :: Vs    ! source voltage

real( c_double ),allocatable :: Vli(:)
real( c_double ),allocatable :: Vlr(:)
real( c_double ),allocatable :: Vri(:)
real( c_double ),allocatable :: Vrr(:)

real( c_double ) :: V,V1,Vnz
  
! START

! set up TLM

dl=0.01d0     ! 1cm mesh
nz=133        ! number of cells=133 i.e. a 1.33m transmission line

Z0=635d0      ! characteristic impedance of the transmission line 
c0=2.0d8      ! velocity of waves on the transmission line 

dt=dl/c0      ! TLM timestep

eps=dt/100d0  ! small time period wrt dt used to determine when to stop the ngspice solution

!nt=2000       ! number of TLM timesteps
nt=2000       ! number of TLM timesteps

tmax=nt*dt    ! maximum time

Zs=100d0       !`source resistance

Vs=1d0        ! source voltage

write(*,*)'dl  =',dl
write(*,*)'nz  =',nz
write(*,*)'len =',nz*dl
write(*,*)'Z0  =',Z0
write(*,*)'c0  =',c0
write(*,*)'dt  =',dt
write(*,*)'nt  =',nt
write(*,*)'tmax=',tmax
write(*,*)'Zs  =',Zs
write(*,*)'Vs  =',Vs

ALLOCATE( Vli(1:nz) )
ALLOCATE( Vlr(1:nz) )
ALLOCATE( Vri(1:nz) )
ALLOCATE( Vrr(1:nz) )

Vli(1:nz)=0d0
Vlr(1:nz)=0d0
Vri(1:nz)=0d0
Vrr(1:nz)=0d0

write(*,*)'Create the Spice_circuit.cir file from the template'

write(Z0_string,'(ES16.6)')Z0
write(dt_string,'(ES16.6)')dt/20d0            ! set the ngspice timestep to something much less than the TLM timestep
write(tmax_string,'(ES16.6)')tmax*1.001d0 +dt ! add an extra bit of time to avoid problems with numerical precision

open(unit=30,file='Spice_circuit_TEMPLATE.cir',status='OLD',ERR=9000)
open(unit=31,file='Spice_circuit.cir')

10 CONTINUE
  read(30,'(A256)',ERR=9010,END=1000)line
  
! look for the string '#Z0_TLM' and replace with Z0
  CALL replace_in_string(line,'#Z0_TLM',Z0_string)
  
! look for the string '#dt_ngspice' and replace with dt
  CALL replace_in_string(line,'#dt_ngspice',dt_string)
  
! look for the string '#tmax_ngspice' and replace with tmax
  CALL replace_in_string(line,'#tmax_ngspice',tmax_string)
  
  write(31,'(A)')trim(line)
  
  GOTO 10  ! read next line

1000 CONTINUE ! jump here when we have reached the end of the input file

close(unit=30)
close(unit=31)

write(*,*)'CALLING ngspice_wrapper_load_circuit'

! The c function needs to be edited to include the following:
! TLM characteristic impedance
! TLM load circuit element(s)
! Transient solution tmax and timestep (.tran) 

istat = ngspice_wrapper_load_circuit( )

write(*,*)'Exit status',istat

n_nodes=0   ! the number of nodes linking TLM and ngspice is not known until the first initdata function call so
            ! set to zero for now.

open(unit=10,file='TLM_solution.dat')

t_F90=0d0  ! initial time

! loop over TLM timestep
do i=1,nt

! SCATTERING PROCESS

! Internal TLM scatter inside cells 
  do iz=1,nz

    V=(2d0*Vli(iz)/Z0+2d0*Vri(iz)/Z0)/(1d0/Z0+1d0/Z0)
    Vlr(iz)=V-Vli(iz)
    Vrr(iz)=V-Vri(iz)
    
  end do
  
! UPDATE THE NGSPICE CIRCUIT SOLUTION
! work out the next time to halt the ngspice solution and write the breakpoint to ngspice

  break_time=dt*dble(i)
  
! Transfer the TLM incident voltage pulse(s) to the ngspice circuit
  do nodeloop=1,n_nodes
  
    opnode=node_list(nodeloop)
    
    if (opnode.EQ.1) then
! source end model
      command_string=''
      write(command_string,'(A,ES16.6)')"alter vtlm1 = ",2d0*Vlr(1)
      istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR);
!      write(*,*)trim(command_string)
      
    else if (opnode.EQ.2) then
! load end model    
      command_string=''
      write(command_string,'(A,ES16.6)')"alter vtlm2 = ",2d0*Vrr(nz)
      istat = ngSpice_Command(trim(command_string)//C_NULL_CHAR);
!      write(*,*)trim(command_string)
      
    end if
    
  end do

! single step through the ngspice solution until the time is greater than or equal to the next break time
  do while((break_time-t_F90).GT.eps) 
  
    istat = ngspice_wrapper_step(t_F90,V_F90,n_nodes,node_list,V_array_F90)
    
  end do
  
! Transfer the ngspice voltage(s) to the TLM solution   
  do nodeloop=1,n_nodes
  
    opnode=node_list(nodeloop)
    
    if (opnode.EQ.1) then
      V1=V_array_F90(nodeloop)
!      write(*,*)'V1 =',V1
    else if (opnode.EQ.2) then
      Vnz=V_array_F90(nodeloop)
!      write(*,*)'Vnz=',Vnz,V_F90
    end if
    
  end do

! CONNECT PROCESS (INCLUDING TERMINATION CONDITIONS)

  time=(i-1)*dt

! TLM connect at source end - this requires the link to ngspice
  
  Vli(1)=V1-Vlr(1)

! Internal TLM connection between cells  
  do iz=1,nz-1

    V=(2d0*Vrr(iz)/Z0+2d0*Vlr(iz+1)/Z0)/(1d0/Z0+1d0/Z0)
    Vri(iz)=V-Vrr(iz)
    Vli(iz+1)=V-Vlr(iz+1)
    
  end do

! TLM connect at load end - this requires the link to ngspice
  
! connect at the load end using ngspice voltage
  Vri(nz)=Vnz-Vrr(nz)

! write voltages at the source and load end to file    
  write(*,*)(i-1)*dt,V1,Vnz
  write(10,*)(i-1)*dt,V1,Vnz

end do

close(unit=10)

DEALLOCATE( Vli )
DEALLOCATE( Vlr )
DEALLOCATE( Vri )
DEALLOCATE( Vrr )

STOP 0

9000 write(*,*)'ERROR: cannot open the template circuit file: Spice_circuit_TEMPLATE.cir'
STOP 1

9010 write(*,*)'ERROR reading template circuit file: Spice_circuit_TEMPLATE.cir'
STOP 1

END SUBROUTINE run_ngspice_TLM
!
! ______________________________________________________________________________
!
!
SUBROUTINE replace_in_string(line,find_string,replace_string)

IMPLICIT NONE

character(LEN=256),INTENT(INOUT) :: line
character(LEN=*),INTENT(IN)      :: find_string
character(LEN=256),INTENT(IN)    :: replace_string

! local variables

integer :: found_pos
integer :: len_line
integer :: len_found
character(LEN=256) :: left_string
character(LEN=256) :: right_string
character(LEN=256) :: opline

! START

!write(*,*)''
!write(*,*)'CALLED:replace_in_string'
!write(*,*)'LINE   :',trim(line)
!write(*,*)'FIND   :',trim(find_string)
!write(*,*)'REPLACE:',trim(replace_string)

len_line=len(trim(line))
len_found=len(trim(find_string))
found_pos=index(line,trim(find_string))

do while (found_pos.NE.0) 

! replace the found string with the replace string
  left_string=''
  left_string=line(1:found_pos-1)
  right_string=''
  right_string=line(found_pos+len_found:len_line)
  
  opline=''
  opline=trim(left_string)//trim(replace_string)//trim(right_string)
  
  line=''
  line=opline
  
  found_pos=index(line,trim(find_string))

end do

END SUBROUTINE replace_in_string
