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
! SUBROUTINE run_TLM_periodic_BC
!
! NAME
!     run_TLM_periodic_BC
!
! DESCRIPTION
!     Main update loop for the TLM solution with periodic boundary conditions applied
!     using the method of Lee and Smith, 
!     "An alternative approach for implementing periodic boundary conditions in the FDTD method using multiple unit cells,"
!     IEEE trans AP vol 54, no2, 2006 pp 698-705
!     
! COMMENTS
!     We note that the solution can go unstable if dielectric or magnetic materials are present
!     This is due to the dispersion of the stub loaded node which allows the fastest error to propagate at twice
!     the speed of light along the coordinate axes of the mesh.
!
!     It may be that some low pass filtering in the process may help to reduce this error. 
!
! HISTORY
!
!     started 4/03/2014 CJS
!
!
SUBROUTINE run_TLM_periodic_BC

USE TLM_general
USE mesh
USE TLM_excitation
USE File_information

USE TLM_output   ! temp for testing

IMPLICIT NONE

! local variables

  character*256	:: ipline
  
  integer	:: op_time_period
  
  integer	:: time_count,start_time_count,last_time_count,last_timestep_time
  integer     	:: time_count_rate
  integer     	:: time_count_max
  
  integer	:: n_timesteps_to_finish
  real*8	:: runtime_per_timestep
  integer	:: time_to_finish,time_to_finish_hrs,time_to_finish_min,time_to_finish_sec
  integer	:: time_to_allocate

! START
  
  CALL write_line('CALLED: run_TLM_periodic_BC',0,output_to_screen_flag)
  
  timestepping_output_to_screen_flag=.FALSE.
!  timestepping_output_to_screen_flag=.TRUE.

  if (rank.eq.0) then
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Mesh'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Nx=',nx 
    write(info_file_unit,*)'Ny=',ny 
    write(info_file_unit,*)'Nz=',nz 
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Total number of cells=',nx*ny*nz
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Number of timesteps=',n_timesteps
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Timestep=',dt
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Number of processes= ',np
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'bicubic_warp_flag= ',bicubic_warp_flag
    write(info_file_unit,*)'frequency_scale_flag= ',frequency_scale_flag
    write(info_file_unit,*)'frequency_scale= ',frequency_scale
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'GGI_TLM solution Started:'
    call write_date_and_time(info_file_unit)
    write(info_file_unit,*)'' 
    write(*,*)'Problem name:',trim(problem_name)
    write(*,*)''
    write(*,*)'TLM solution Started:'
    call write_date_and_time(0)
    write(*,*)'' 
  end if

#if defined(MPI)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
#endif

  if (rank.eq.0) then
  
#if defined(SEQ)
    call system("ps u -C GGI_TLM_SEQ > GGI_TLM_memory_usage.txt ")
#elif defined(MPI)
    call system("ps u -C GGI_TLM_MPI > GGI_TLM_memory_usage.txt ")
#endif
    
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,'(A)')""
    write(info_file_unit,'(A)')"Memory Usage:"
    write(info_file_unit,'(A)')""

    write(*,'(A)')'____________________________________________________'
    write(*,'(A)')""
    write(*,'(A)')"Memory Usage:"
    write(*,'(A)')""
    
    open(unit=scratch_file_unit,file='GGI_TLM_memory_usage.txt')
5   CONTINUE
    read(scratch_file_unit,'(A256)',end=6)ipline
    write(info_file_unit,'(A)')trim(ipline)
    write(*,'(A)')trim(ipline)
    
    GOTO 5
    
6   CONTINUE

    write(info_file_unit,'(A)')""
    write(info_file_unit,'(A)')'____________________________________________________'
    write(info_file_unit,'(A)')""
    write(*,'(A)')""
    write(*,'(A)')'____________________________________________________'
    write(*,'(A)')""
    
  end if
    
  CALL system_clock(start_time_count,time_count_rate,time_count_max)
  last_time_count=start_time_count
  last_timestep_time=1

  do timestep=1,n_timesteps

    op_time_period=min(100,10**INT(log10(dble(timestep))))

    if (rank.eq.0) then 

      if (timestep.EQ.1) then
      
	write(6,8000,advance='no')'Timestep ',timestep,' of ',n_timesteps
	flush(6)
	
      else if (timestep.EQ.n_timesteps) then
            
	write(6,'(A)',advance='no')char(13)
	write(6,8000)'Timestep ',timestep,' of ',n_timesteps
	flush(6)    
	  
      else if ( mod(timestep,op_time_period).EQ.0 ) then
      
! estimate time to finish
      
        CALL system_clock(time_count,time_count_rate,time_count_max)
    
        runtime_per_timestep=dble(time_count-last_time_count)/dble((timestep-last_timestep_time)*time_count_rate)
        
        n_timesteps_to_finish=n_timesteps+1-timestep
        time_to_finish=NINT( dble(n_timesteps_to_finish)*runtime_per_timestep )
			      
        time_to_finish_hrs=INT(time_to_finish/3600)
        time_to_allocate=time_to_finish-time_to_finish_hrs*3600
        time_to_finish_min=INT(time_to_allocate/60)
        time_to_finish_sec=time_to_allocate-time_to_finish_min*60
      
	write(6,'(A)',advance='no')char(13)
	write(6,8010,advance='no')'Timestep ',timestep,' of ',n_timesteps,	&
	             '  Estimated time to finish: ',time_to_finish_hrs,':',time_to_finish_min,':',time_to_finish_sec
	flush(6)
    
8000    format(A9,I10,A4,I10)    
8010    format(A9,I10,A4,I10,A27,I6.2,A,I2.2,A,I2.2)    
	  
      end if
      
    end if ! rank=0

    time=(timestep-1)*dt
    
    CALL excitation()
  
    CALL scatter()
  
    CALL cell_output()
    
    time=time+dt/2d0
    
    if (np.gt.1) then 
    
      CALL TLM_parallel_pass_data_1()
      
      CALL Cable_parallel_pass_data_1()
      
    end if
    
    CALL mode_stir_surfaces()

    CALL connect()
  
    CALL face_output()
    
    CALL cable_output()
    
    CALL outer_boundary()
    
    CALL wrap_outer_boundary()
    
    if (np.gt.1) then 
    
      CALL TLM_parallel_pass_data_2()
      
      CALL Cable_parallel_pass_data_2()
      
    end if
  
  end do ! next timestep

  if (rank.eq.0) then
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'TLM solution Finished:'
    call write_date_and_time(info_file_unit)
    write(info_file_unit,*)'' 
    
    write(*,*)'' 
    write(*,*)'TLM solution Finished:'
    call write_date_and_time(0)
    write(*,*)'' 
    
! calculate runtime
      
    CALL system_clock(time_count,time_count_rate,time_count_max)
    
    runtime_per_timestep=dble(time_count-start_time_count)/dble(n_timesteps*time_count_rate)

    time_to_allocate=(time_count-start_time_count)/time_count_rate
			      
    time_to_finish_hrs=INT(time_to_allocate/3600)
    time_to_allocate=time_to_allocate-time_to_finish_hrs*3600
    time_to_finish_min=INT(time_to_allocate/60)
    time_to_finish_sec=time_to_allocate-time_to_finish_min*60
    
    write(info_file_unit,8020)'Run time: ',time_to_finish_hrs,':',time_to_finish_min,':',time_to_finish_sec
    write(info_file_unit,*)'' 
    write(info_file_unit,8030)'Run time per timestep: ',runtime_per_timestep,' seconds'
    write(info_file_unit,*)'' 
    
    write(*,8020)'Run time: ',time_to_finish_hrs,':',time_to_finish_min,':',time_to_finish_sec
    write(*,*)'' 
    write(*,8030)'Run time per timestep: ',runtime_per_timestep,' seconds'
    write(*,*)'' 
    
8020  format(A10,I6.2,A,I2.2,A,I2.2)    
8030  format(A23,E10.2,A8)    
    
  end if
  
  CALL write_line('FINISHED: run_TLM_periodic_BC',0,output_to_screen_flag)

  RETURN

END SUBROUTINE run_TLM_periodic_BC
!
! SUBROUTINE setup_periodic_bc
!
! Name setup_periodic_bc
!     
!
! Description
!     work out the timestep loop for periodic structures
!
! Comments:
!      
!
! History
!
!     started 7/10/10 CJS
!     

SUBROUTINE setup_periodic_bc()

USE mesh
USE constants
USE TLM_general
USE TLM_periodic

USE cell_parameters
USE file_information
USE TLM_excitation

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  real*8 lxbc,lybc
  real*8 ld,lp
  real*8 theta_w,theta_g,theta_d
  real*8 kx,ky,kz
  real*8 tdx,tdy,td
  real*8 tau_d,taud_x,taud_y,tau_e
  integer ntd,nterr
  real*8 t_0,terr_x,terr_y,terr
  real*8 tmax,tmin
  real*8 r,phi,theta
  
  integer nt1,nt2

! START

  write(*,*)'CALLED: setup_periodic_bc'

! stage 1. Look for a single huygens surface and get k vector

  if (n_huygens_surfaces.ne.1) then
    write(*,*)'Error: we need a single Huygens surface for periodic structure process'
    stop
  end if
  
! stage 2. calculate delays in x and y directions

! calculate incident vector
  r=1d0
  theta=huygens_surface%Ktheta
  phi=huygens_surface%Kphi

!!!! OLD NOTATION    
!  call rthetaphi_to_xyz(r,theta,phi,kx,ky,kz )
  write(info_file_unit,*)'Huygens surface, theta=',theta,theta*180d0/pi
  write(info_file_unit,*)'Huygens surface, phi  =',phi,phi*180d0/pi
    
! tweak the angle so that delays are an integer number of timesteps
! in the periodic structure update  

  lxbc=(nx/2)*dl
  lybc=(ny/2)*dl
  
  tdx=kx*lxbc/c0
  tdy=ky*lybc/c0
  td=tdx+tdy

  write(info_file_unit,*)'lxbc=',lxbc
  write(info_file_unit,*)'lybc=',lybc

  write(info_file_unit,*)'tdx=',tdx
  write(info_file_unit,*)'tdy=',tdy
  
  ntdx=Nint(tdx/dt)
  ntdy=Nint(tdy/dt)
  ntd=ntdx+ntdy
  
  write(info_file_unit,*)'ntdx=',ntdx
  write(info_file_unit,*)'ntdy=',ntdy
  
! reset Huygens surface data
  tdx=ntdx*dt
  tdy=ntdy*dt
  
  kx=tdx*c0/lxbc
  ky=tdy*c0/lybc
  kz=sqrt(1d0-kx*kx-ky*ky)

! OLD NOTATION  
!!!!  call xyz_to_rthetaphi( kx,ky,kz,r,theta,phi )
  huygens_surface%Ktheta=theta
  huygens_surface%Kphi=phi
  
  write(info_file_unit,*)'Revised Huygens surface, theta=',theta,theta*180d0/pi
  write(info_file_unit,*)'Revised Huygens surface, phi  =',phi,phi*180d0/pi
  
  if (kx.ge.0d0) then
    pbc_xface=face_xmax
  else 
    pbc_xface=face_xmin
  end if
  
  if (ky.ge.0d0) then
    pbc_yface=face_ymax
  else 
    pbc_yface=face_ymin
  end if
  
  write(info_file_unit,*)'kx=',kx
  write(info_file_unit,*)'ky=',ky
  write(info_file_unit,*)'kz=',kz
  
  write(info_file_unit,*)'pbc_xface=',pbc_xface
  write(info_file_unit,*)'pbc_yface=',pbc_yface
  
  if (pbc_xface.eq.face_xmin) then
    write(*,*)' Require kx.gt.0 for PBC'
    write(*,*)'kx=',kx
    write(*,*)'ky=',ky
    write(*,*)'kz=',kz
    stop
  end if
  if (pbc_yface.eq.face_ymin) then
    write(*,*)' Require ky.gt.0 for PBC'
    write(*,*)'kx=',kx
    write(*,*)'ky=',ky
    write(*,*)'kz=',kz
    stop
  end if
     
  write(info_file_unit,*)'ntdx=',ntdx
  write(info_file_unit,*)'ntdy=',ntdy
  write(info_file_unit,*)'ntd =',ntdx+ntdy
  
  theta_w=atan2(ky,kx)
  theta_g=atan2(lybc,lxbc)
  theta_d=theta_w-theta_g
  
  ld=sqrt(lxbc*lxbc+lybc*lybc)
  lp=ld*cos(theta_d)

  t_0=(kx*lxbc+ky*lybc)/c0
  
  write(*,*)'t_0',t_0
 
  terr_x=(lxbc-dl*0d0)/c0 
  terr_y=(lybc-dl*0d0)/c0
  
  terr=min(terr_x,terr_y)

! old  
!  nterr=int(terr/dt)
  
! new
  nterr=int(terr/dt)-3    ! can remove another timestep - botched could maybe change back if we copy data
                          ! after scatter and before connect stage
  terr=nterr*dt
  
  write(info_file_unit,*)'terr_x',terr_x
  write(info_file_unit,*)'terr_y',terr_y
  write(info_file_unit,*)'terr  ',terr
  write(info_file_unit,*)'nterr ',nterr
  
  taud_x=kx*lxbc/c0
  
  tau_e=terr
  tau_d=td/c0
  
  write(info_file_unit,*)'Inequality'
  write(info_file_unit,*)lxbc*kx+lybc*ky,' < ',min(lxbc,lybc)
  
  if (tau_d.ge.tau_e) then
    write(info_file_unit,*)'tau_d.ge.tau_e'
    write(info_file_unit,*)'tau_e=',tau_e
    write(info_file_unit,*)'tau_d=',tau_d
    stop
  end if

  write(*,*)
  
  tmin=0.0
  tmax=simulation_time
  write(info_file_unit,*)'Simulation time ',tmax
 
! stage 3. Calculate the run loop parameters
  
  t_cycle=terr-tdx-tdy
  tmax_cycle=terr
  nt_cycle=nint(t_cycle/dt)
  ntmax_cycle=nint(tmax_cycle/dt)
  
  write(info_file_unit,*)'t_cycle=',t_cycle,' nt_cycle=',nt_cycle,ntmax_cycle-ntdx-ntdy
  write(info_file_unit,*)'tmax_cycle=',tmax_cycle,' ntmax_cycle=',ntmax_cycle
  
  n_cycles=int((tmax+t_0)/t_cycle)+1
  write(info_file_unit,*)'n_cycles=',n_cycles
  
  ntmin=nint(-t_0/dt) ! time for wave to reach region 1
  ntmax=nint(tmax/dt)
  
  nt_run(1)=ntmax_cycle-ntdx-ntdy
  n_save(1)=1
  
  if (ntdx.gt.ntdy) then
    nt_run(2)=ntmax_cycle-ntdx
    n_save(2)=2
    nt_run(3)=ntmax_cycle-ntdy
    n_save(3)=3
  else
    nt_run(2)=ntmax_cycle-ntdy
    n_save(2)=3
    nt_run(3)=ntmax_cycle-ntdx
    n_save(3)=2
  end if
  nt_run(4)=ntmax_cycle
  n_save(4)=4
  
  n_pbc_timesteps(1)=nt_run(1)
  n_pbc_timesteps(2)=nt_run(2)-nt_run(1)
  n_pbc_timesteps(3)=nt_run(3)-nt_run(1)
  n_pbc_timesteps(4)=nt_run(4)-nt_run(1)
  
  
  write(info_file_unit,*)'t_cycle =',t_cycle
  write(info_file_unit,*)'nt_cycle =',nt_cycle
  write(info_file_unit,*)'ntmin=',ntmin
  write(info_file_unit,*)'ntmax=',ntmax
  write(info_file_unit,*)'nt_run1 =',nt_run(1),' n_save(1)',n_save(1)
  write(info_file_unit,*)'nt_run2 =',nt_run(2),' n_save(2)',n_save(2)
  write(info_file_unit,*)'nt_run3 =',nt_run(3),' n_save(3)',n_save(3)
  write(info_file_unit,*)'nt_run4 =',nt_run(4),' n_save(4)',n_save(4)
  
  nt_cycle_start=ntmin
 
  write(*,*)'FINISHED: setup_periodic_bc'
  
  return
  
END SUBROUTINE setup_periodic_bc
