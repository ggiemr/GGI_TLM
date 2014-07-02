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
USE TLM_periodic
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
  
  real*8	:: t1,t2
  integer	:: run,cycle,count,iteration,n
  integer	:: nt1,nt2

! START
  
  CALL write_line('CALLED: run_TLM_periodic_BC',0,output_to_screen_flag)
  
  timestepping_output_to_screen_flag=.FALSE.
!  timestepping_output_to_screen_flag=.TRUE.

! This call is now made from GGI_TLM.F90 before the Huygens surface is set up
!  CALL setup_periodic_bc()
  
  psox=nx/2
  psoy=ny/2
  psoz=0

! Allocate memory for periodic boundary data

  allocate( V_pbc_save(1:12,1:nx,1:ny,nzmin:nzmax) )
  
  n_xbc_point=ntdx*2
  n_ybc_point=ntdy*2
  
  allocate( V_pbc_x(1:ny,nz1:nz2,1:2,1:n_xbc_point) )
  allocate( V_pbc_y(1:nx,nz1:nz2,1:2,1:n_ybc_point) )
  allocate( V_pbc_x_save(1:ny,nz1:nz2,1:2,1:n_xbc_point) )
  allocate( V_pbc_y_save(1:nx,nz1:nz2,1:2,1:n_ybc_point) )
  
  V_pbc_save(1:12,1:nx,1:ny,nz1:nz2)=0d0
  V_pbc_x(1:ny,nz1:nz2,1:2,1:n_xbc_point)=0d0
  V_pbc_y(1:nx,nz1:nz2,1:2,1:n_ybc_point)=0d0
  V_pbc_x_save(1:ny,nz1:nz2,1:2,1:n_xbc_point)=0d0
  V_pbc_y_save(1:nx,nz1:nz2,1:2,1:n_ybc_point)=0d0

  xbc_write_point_ymin=1
  xbc_write_point_ymax=1
  ybc_write_point_xmin=1
  ybc_write_point_xmax=1

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

! Timestep loop(s)
  
  count=0
  iteration=0
!  timestep=0
  
  write(6,8000,advance='no')'Timestep ',timestep,' of ',n_timesteps
8000    format(A9,I10,A4,I10)    
  flush(6)
  
  do cycle=1,n_cycles

    do run=1,4
    
! if run=1 then this part of the simulation gives accurate results in region 1
! and results should be recorded if t>=0

      if (run.eq.1) then
      
        nt1=nt_cycle_start
	call TLM_pbc_copy_saved_field()
	
      else
      
        nt1=nt_cycle_start+nt_run(run-1)
	
      end if
      
      nt2=nt_cycle_start+nt_run(run)
   
      t1=nt1*dt
      t2=nt2*dt    
      
      do n=nt1+1,nt2
      
        timestep=n
        time=dt*(n-1)
	
        count=count+1
  
        if (run.eq.1) then	
! UPDATE WITH OUTPUT
      
          if (rank.eq.0) then
	    write(6,'(A)',advance='no')char(13)
	    write(6,8010,advance='no')'Timestep ',timestep,' of ',n_timesteps,' : stage=',run,' of 4 : cycle=',cycle,' of ',n_cycles
	    flush(6)   
8010        format(A9,I10,A4,I10,A9,I2,A14,I6,A4,I6)    
	  end if
	  
	  if (time.ge.0D0) then
	  
            iteration=iteration+1

	  end if
	      
          CALL excitation()  
          CALL scatter()  
          CALL cell_output()
    
          time=time+dt/2d0
    
          if (np.gt.1) then 
    
            CALL TLM_parallel_pass_data_1()      
            CALL Cable_parallel_pass_data_1()
      
          end if
    
          CALL connect()  
          CALL face_output()   
	  CALL TLM_pbc_outer_boundary
    
          if (np.gt.1) then 
    
            CALL TLM_parallel_pass_data_2()      
            CALL Cable_parallel_pass_data_2()
      
          end if

! END OF UPDATE SEQUENCE WITH OUTPUT
	  
	  if (time+dt/2d0.ge.simulation_time) then ! Finish as the final timestep may be mid cycle. 
	    goto 1000
	  end if
  
 	  	  
        else ! run.ne.1
	
          if (nt1.ne.nt2) then
! only update if we have timesteps here otherwise there is no point
      
            if (rank.eq.0) then
	      write(6,'(A)',advance='no')char(13)
	      write(6,8010,advance='no')'Timestep ',timestep,' of ',n_timesteps,' : stage=',run,' of 4 : cycle=',cycle,' of ',n_cycles
	      flush(6)   
	    end if

! UPDATE WITH NO OUTPUT
    
            CALL excitation()  
            CALL scatter()
    
            time=time+dt/2d0
    
            if (np.gt.1) then 
    
              CALL TLM_parallel_pass_data_1()      
              CALL Cable_parallel_pass_data_1()
      
            end if

            CALL connect()      
	    CALL TLM_pbc_outer_boundary
    
            if (np.gt.1) then 
    
              CALL TLM_parallel_pass_data_2()      
              CALL Cable_parallel_pass_data_2()
      
            end if

! END OF UPDATE SEQUENCE 
	    
	  end if ! (nt1.ne.nt2)
	  	  
        end if ! run.ne.1
	
      end do ! next timestep in this run
  
! Save field in region n_save(run)
      call TLM_pbc_save_field(n_save(run))  
       
      if (run.eq.4) then ! save next start timestep
        nt_cycle_start=nt_cycle_start+nt_run(1)
      end if
      
    end do ! next run 
      
  end do ! next cycle

1000 CONTINUE    ! jump here at last timestep
              
  if (rank.eq.0) then
  
    write(6,'(A)',advance='no')char(13)
    write(6,8000)'Timestep ',timestep,' of ',n_timesteps
    flush(6)    
  
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
!      We may need to check that the polarisation is correct for the case when theta=0 as the process can
!      then set phi=0 whatever it is set to initially.
!
! History
!
!     started 29/4/14 CJS based on the original HIRF-SE code Fieldsolve
!     

SUBROUTINE setup_periodic_bc()

USE mesh
USE geometry_types
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
  type(xyz)	:: k
  real*8 tdx,tdy,td
  real*8 tau_d,taud_x,taud_y,tau_e
  integer ntd,nterr
  real*8 t_0,terr_x,terr_y,terr
  real*8 tmax,tmin
  real*8 r,phi,theta,phi_in
  
  real*8 dxdt,dxdp
  real*8 dydt,dydp
  real*8 dzdt,dzdp
  real*8 vx,vy,vz
  real*8 kx,ky,kz
  real*8 norm
  real*8 Ei(3),Hi(3)
  type(xyz)		:: Kvector
  
  integer nt1,nt2

! START

  write(*,*)'CALLED: setup_periodic_bc'
  
  write(info_file_unit,*)'_______________________________________________'
  write(info_file_unit,*)''
  write(info_file_unit,*)'Periodic Boundary Condition Data'
  write(info_file_unit,*)''

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

! save initial phi so that polarisation ambiguity may be resolved for theta=0  
  phi_in=phi

  CALL rthetaphi_to_xyz_point(r,theta,phi,k)
  
  write(info_file_unit,*)'Huygens surface, theta=',theta,theta*180d0/pi
  write(info_file_unit,*)'Huygens surface, phi  =',phi,phi*180d0/pi
    
! tweak the angle so that delays are an integer number of timesteps
! in the periodic structure update  

  lxbc=(nx/2)*dl
  lybc=(ny/2)*dl
  
  tdx=k%x*lxbc/c0
  tdy=k%y*lybc/c0
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
  
  k%x=tdx*c0/lxbc
  k%y=tdy*c0/lybc
  k%z=sqrt(1d0-k%x*k%x-k%y*k%y)

  CALL xyz_point_to_rthetaphi(k,r,theta,phi)

! resolve polarisation ambiguity if theta=0  
  if (theta.eq.0d0) then
    phi=phi_in
  end if
  
  huygens_surface%Ktheta=theta
  huygens_surface%Kphi=phi
  
  write(info_file_unit,*)'Revised Huygens surface, theta=',theta,theta*180d0/pi
  write(info_file_unit,*)'Revised Huygens surface, phi  =',phi,phi*180d0/pi
  
! re-calculate other Huygens surface data

! calculate incident vector
  r=1d0
  theta=huygens_surface%Ktheta
  phi=huygens_surface%Kphi
  
  CALL rthetaphi_to_xyz_point(r,theta,phi,Kvector)
  
  kx=Kvector%x
  ky=Kvector%y
  kz=Kvector%z
  
  huygens_surface%Ki(1)=kx
  huygens_surface%Ki(2)=ky
  huygens_surface%Ki(3)=kz
  
  if (rank.eq.0) write(info_file_unit,8000)'K vector=',huygens_surface%Ki(1),huygens_surface%Ki(2),huygens_surface%Ki(3)
8000 format(A,3F12.6)

  dxdt= r*cos(theta)*cos(phi)
  dxdp=-r*sin(theta)*sin(phi)
  dydt= r*cos(theta)*sin(phi)
  dydp= r*sin(theta)*cos(phi)
  dzdt=-r*sin(theta)
  dzdp=0d0
  vx=huygens_surface%Etheta*dxdt+huygens_surface%Ephi*dxdp
  vy=huygens_surface%Etheta*dydt+huygens_surface%Ephi*dydp
  vz=huygens_surface%Etheta*dzdt+huygens_surface%Ephi*dzdp
  
  norm=sqrt(vx*vx+vy*vy+vz*vz)
  
  if (norm.eq.0d0) then
    write(*,*)'Error calculating Huygens surface data'
    write(*,8000)'E vector(theta)=',dxdt,dydt,dzdt
    write(*,8000)'E vector(phi)=',dxdp,dydp,dzdp
    write(*,8000)'E vector=',huygens_surface%Etheta,huygens_surface%Ephi
    STOP
  end if

  Ei(1)=vx/norm
  Ei(2)=vy/norm
  Ei(3)=vz/norm    
!note K cross E = H

  call vector_product(kx,ky,kz,Ei(1),Ei(2),Ei(3),Hi(1),Hi(2),Hi(3))
  Hi(:)=Hi(:)/Z0
  
  huygens_surface%Ei(1)=Ei(1)
  huygens_surface%Ei(2)=Ei(2)
  huygens_surface%Ei(3)=Ei(3)
  huygens_surface%Hi(1)=Hi(1)
  huygens_surface%Hi(2)=Hi(2)
  huygens_surface%Hi(3)=Hi(3)
  
  if (rank.eq.0) write(info_file_unit,8000)'E vector=',huygens_surface%Ei(1),huygens_surface%Ei(2),huygens_surface%Ei(3)
  if (rank.eq.0) write(info_file_unit,8000)'H vector=',huygens_surface%Hi(1),huygens_surface%Hi(2),huygens_surface%Hi(3)
  
  if (k%x.ge.0d0) then
    pbc_xface=face_xmax
  else 
    pbc_xface=face_xmin
  end if
  
  if (k%y.ge.0d0) then
    pbc_yface=face_ymax
  else 
    pbc_yface=face_ymin
  end if
  
  write(info_file_unit,*)'kx=',k%x
  write(info_file_unit,*)'ky=',k%y
  write(info_file_unit,*)'kz=',k%z
  
  write(info_file_unit,*)'pbc_xface=',pbc_xface
  write(info_file_unit,*)'pbc_yface=',pbc_yface
  
  if (pbc_xface.eq.face_xmin) then
    write(*,*)' Require kx.gt.0 for PBC'
    write(*,*)'kx=',k%x
    write(*,*)'ky=',k%y
    write(*,*)'kz=',k%z
    stop
  end if
  if (pbc_yface.eq.face_ymin) then
    write(*,*)' Require ky.gt.0 for PBC'
    write(*,*)'kx=',k%x
    write(*,*)'ky=',k%y
    write(*,*)'kz=',k%z
    stop
  end if
     
  write(info_file_unit,*)'ntdx=',ntdx
  write(info_file_unit,*)'ntdy=',ntdy
  write(info_file_unit,*)'ntd =',ntdx+ntdy
  
  theta_w=atan2(k%y,k%x)
  theta_g=atan2(lybc,lxbc)
  theta_d=theta_w-theta_g
  
  ld=sqrt(lxbc*lxbc+lybc*lybc)
  lp=ld*cos(theta_d)

  t_0=(k%x*lxbc+k%y*lybc)/c0
  
  write(info_file_unit,*)'t_0',t_0
 
  terr_x=(lxbc-dl*0d0)/c0 
  terr_y=(lybc-dl*0d0)/c0
  
  terr=min(terr_x,terr_y)
  
! new
  nterr=int(terr/dt)-3    ! can remove another timestep - botched could maybe change back if we copy data
                          ! after scatter and before connect stage
  terr=nterr*dt
  
  write(info_file_unit,*)'terr_x',terr_x
  write(info_file_unit,*)'terr_y',terr_y
  write(info_file_unit,*)'terr  ',terr
  write(info_file_unit,*)'nterr ',nterr
  
  taud_x=k%x*lxbc/c0
  
  tau_e=terr
  tau_d=td/c0
  
  write(info_file_unit,*)'Inequality'
  write(info_file_unit,*)lxbc*k%x+lybc*k%y,' < ',min(lxbc,lybc)
  
  if (tau_d.ge.tau_e) then
    write(info_file_unit,*)'tau_d.ge.tau_e'
    write(info_file_unit,*)'tau_e=',tau_e
    write(info_file_unit,*)'tau_d=',tau_d
    stop
  end if

  write(info_file_unit,*)
  
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

!
! Name TLM_pbc_outer_boundary
!     
!
! Description
!     
!
! Comments:
!      boundary conditions for periodic boundary implementation 
!
! History
!
!     started 29/4/14 CJS based on the original HIRF-SE code Fieldsolve
!     
!

SUBROUTINE TLM_pbc_outer_boundary()

USE mesh
USE constants
USE TLM_general
USE TLM_periodic
USE cell_parameters

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer ix,iy,iz
  
  integer :: xbc_read_point_ymin
  integer :: ybc_read_point_xmin
  
  integer :: xbc_read_point_ymax
  integer :: ybc_read_point_xmax

! START

!  write(*,*)'CALLED: outer_boundary'      

!  write(temp_file_unit,*)'Called : outer boundary'

! faces normal to x i.e. the y z plane
    do iy=1,ny
      do iz=nz1,nz2
    
          ix=1
          V(Vy_xmin,ix,iy,iz)=0d0 
          V(Vz_xmin,ix,iy,iz)=0d0

! no boundary condition on xmax	
	
      end do  ! next z cell
    end do    ! next y cell
  
! faces normal to y i.e. the x z plane
    do ix=1,nx
      do iz=nz1,nz2
    
          iy=1
          V(Vx_ymin,ix,iy,iz)=0d0
          V(Vz_ymin,ix,iy,iz)=0d0
	
! no boundary condition on ymax	
	  
      end do  ! next z cell
    end do    ! next x cell
    
! zmin face
  iz=1
  
  if (rank.eq.0) then
  
    do iy=1,ny
      do ix=1,nx
      
        V(Vx_zmin,ix,iy,iz)=V(Vx_zmin,ix,iy,iz)*R_zmin
        V(Vy_zmin,ix,iy,iz)=V(Vy_zmin,ix,iy,iz)*R_zmin
  	  
      end do  ! next x cell
    end do    ! next y cell
    
  end if  
    
! zmax face
  iz=nz
  
  if (rank.eq.np-1) then
    
    do iy=1,ny
      do ix=1,nx
      
        V(Vx_zmax,ix,iy,iz)=V(Vx_zmax,ix,iy,iz)*R_zmax
        V(Vy_zmax,ix,iy,iz)=V(Vy_zmax,ix,iy,iz)*R_zmax
  	  
      end do  ! next x cell
    end do    ! next y cell
    
  end if  

          
! PBC stuff

  xbc_write_point_ymin=xbc_write_point_ymin+1
  if (xbc_write_point_ymin.gt.n_xbc_point) xbc_write_point_ymin=1
  
  xbc_read_point_ymin=xbc_write_point_ymin-ntdx
  if (xbc_read_point_ymin.le.0) xbc_read_point_ymin=xbc_read_point_ymin+n_xbc_point
  
  xbc_write_point_ymax=xbc_write_point_ymax+1
  if (xbc_write_point_ymax.gt.n_xbc_point) xbc_write_point_ymax=1
  
  xbc_read_point_ymax=xbc_write_point_ymax-ntdx
  if (xbc_read_point_ymax.le.0) xbc_read_point_ymax=xbc_read_point_ymax+n_xbc_point
  
!  write(*,*)'X Write points',xbc_write_point_ymin,xbc_write_point_ymax
  
  if (pbc_xface.eq.face_xmax) then
  
    if (ntdx.ne.0) then      
    
      V(Vy_xmax,nx,1:ny/2,   nz1:nz2)=V_pbc_x(1:ny/2,   nz1:nz2,1,xbc_read_point_ymin)
      V(Vy_xmax,nx,ny/2+1:ny,nz1:nz2)=V_pbc_x(ny/2+1:ny,nz1:nz2,1,xbc_read_point_ymax)  

      V(Vz_xmax,nx,1:ny/2,   nz1:nz2)=V_pbc_x(1:ny/2   ,nz1:nz2,2,xbc_read_point_ymin)      
      V(Vz_xmax,nx,ny/2+1:ny,nz1:nz2)=V_pbc_x(ny/2+1:ny,nz1:nz2,2,xbc_read_point_ymax)    

      V_pbc_x(1:ny/2,   nz1:nz2,1,xbc_write_point_ymin)=V(Vy_xmax,nx/2,1:ny/2,    nz1:nz2)
      V_pbc_x(ny/2+1:ny,nz1:nz2,1,xbc_write_point_ymax)=V(Vy_xmax,nx/2,ny/2+1:ny, nz1:nz2)
      
      V_pbc_x(1:ny/2,   nz1:nz2,2,xbc_write_point_ymin)=V(Vz_xmax,nx/2,1:ny/2,    nz1:nz2)
      V_pbc_x(ny/2+1:ny,nz1:nz2,2,xbc_write_point_ymax)=V(Vz_xmax,nx/2,ny/2+1:ny, nz1:nz2)
      
    else
    
      V(Vy_xmax,nx,1:ny,nz1:nz2)=V(Vy_xmax,nx/2,1:ny,nz1:nz2)
      V(Vz_xmax,nx,1:ny,nz1:nz2)=V(Vz_xmax,nx/2,1:ny,nz1:nz2)
      
    end if 
    
  else if (pbc_xface.eq.face_xmin) then
  
    write(*,*)'Error in TLM_pbc_outer_boundary: pbc_xface.eq.face_xmin'
    stop
   
  end if 

  ybc_write_point_xmin=ybc_write_point_xmin+1
  if (ybc_write_point_xmin.gt.n_ybc_point) ybc_write_point_xmin=1
  
  ybc_read_point_xmin=ybc_write_point_xmin-ntdy
  if (ybc_read_point_xmin.le.0) ybc_read_point_xmin=ybc_read_point_xmin+n_ybc_point
  
  ybc_write_point_xmax=ybc_write_point_xmax+1
  if (ybc_write_point_xmax.gt.n_ybc_point) ybc_write_point_xmax=1
  
  ybc_read_point_xmax=ybc_write_point_xmax-ntdy
  if (ybc_read_point_xmax.le.0) ybc_read_point_xmax=ybc_read_point_xmax+n_ybc_point
  
  if (pbc_yface.eq.face_ymax) then
  
    if (ntdy.ne.0) then
    
      V(Vx_ymax,1:nx/2,   ny,nz1:nz2)=V_pbc_y(1:nx/2,   nz1:nz2,1,ybc_read_point_xmin)
      V(Vx_ymax,nx/2+1:nx,ny,nz1:nz2)=V_pbc_y(nx/2+1:nx,nz1:nz2,1,ybc_read_point_xmax)
      
      V(Vz_ymax,1:nx/2,   ny,nz1:nz2)=V_pbc_y(1:nx/2,   nz1:nz2,2,ybc_read_point_xmin)
      V(Vz_ymax,nx/2+1:nx,ny,nz1:nz2)=V_pbc_y(nx/2+1:nx,nz1:nz2,2,ybc_read_point_xmax)
            
      V_pbc_y(1:nx/2,   nz1:nz2,1,ybc_write_point_xmin)=V(Vx_ymax,1:nx/2,   ny/2,nz1:nz2)
      V_pbc_y(nx/2+1:nx,nz1:nz2,1,ybc_write_point_xmax)=V(Vx_ymax,nx/2+1:nx,ny/2,nz1:nz2)

      V_pbc_y(1:nx/2,   nz1:nz2,2,ybc_write_point_xmin)=V(Vz_ymax,1:nx/2,   ny/2,nz1:nz2)
      V_pbc_y(nx/2+1:nx,nz1:nz2,2,ybc_write_point_xmax)=V(Vz_ymax,nx/2+1:nx,ny/2,nz1:nz2)
            
    else
    
      V(Vx_ymax,1:nx,ny,nz1:nz2)=V(Vx_ymax,1:nx,ny/2,nz1:nz2)
      V(Vz_ymax,1:nx,ny,nz1:nz2)=V(Vz_ymax,1:nx,ny/2,nz1:nz2)
      
    end if
    
  else if (pbc_yface.eq.face_ymin) then
     
    write(*,*)'Error in TLM_pbc_outer_boundary: pbc_xface.eq.face_ymin'
    stop
 
  end if 

!  write(*,*)'FINISHED: TLM_pbc_outer_boundary'
  return
  
END SUBROUTINE TLM_pbc_outer_boundary
!
! SUBROUTINE TLM_pbc_save_field
!
! Name TLM_pbc_save_field
!     
!
! Description
!     copy mesh data for periodic structures
!
! Comments:
!     We need to include material information here too
!
! History
!
!     started 29/4/14 CJS based on the original HIRF-SE code Fieldsolve
!     

SUBROUTINE TLM_pbc_save_field(region)

USE mesh
USE constants
USE TLM_general
USE TLM_periodic
USE cell_parameters
USE TLM_volume_materials
USE TLM_surface_materials

IMPLICIT NONE

! variables passed to subroutine

  integer :: region
  integer :: cx,cy,cz,face
  integer :: point,point2,mat,i

! local_variables

! START

!  write(*,*)'CALLED: TLM_pbc_save_field'

! copy volume material data  
  
  if (region.eq.1) then
  
    V_pbc_save(1:12,nx/2+1:nx,ny/2+1:ny,nz1:nz2)=V(1:12,nx/2+1:nx,ny/2+1:ny,nz1:nz2)

!save boundary stuff here
    V_pbc_x_save(ny/2+1:ny,nz1:nz2,1:2,1:n_xbc_point)=V_pbc_x(ny/2+1:ny,nz1:nz2,1:2,1:n_xbc_point)
    xbc_write_point_save_ymax=xbc_write_point_ymax

    V_pbc_y_save(nx/2+1:nx,nz1:nz2,1:2,1:n_ybc_point)=V_pbc_y(nx/2+1:nx,nz1:nz2,1:2,1:n_ybc_point)
    ybc_write_point_save_xmax=ybc_write_point_xmax
    
  else if (region.eq.2) then
  
    V_pbc_save(1:12,nx/2+1:nx,1:ny/2,nz1:nz2)=V(1:12,nx/2+1:nx,ny/2+1:ny,nz1:nz2)

!save boundary stuff here
    V_pbc_x_save(1:ny/2,nz1:nz2,1:2,1:n_xbc_point)=V_pbc_x(ny/2+1:ny,nz1:nz2,1:2,1:n_xbc_point)
    xbc_write_point_save_ymin=xbc_write_point_ymax
    
  else if (region.eq.3) then
  
    V_pbc_save(1:12,1:nx/2,ny/2+1:ny,nz1:nz2)=V(1:12,nx/2+1:nx,ny/2+1:ny,nz1:nz2)

!save boundary stuff here
    V_pbc_y_save(1:nx/2,nz1:nz2,1:2,1:n_ybc_point)=V_pbc_y(nx/2+1:nx,nz1:nz2,1:2,1:n_ybc_point)
    ybc_write_point_save_xmin=ybc_write_point_xmax
    
  else if (region.eq.4) then
  
    V_pbc_save(1:12,1:nx/2,1:ny/2,nz1:nz2)=V(1:12,nx/2+1:nx,ny/2+1:ny,nz1:nz2)
  
  end if

! need to save dispersive material stuff here
	  
! need to save volume material stuff here
  
!  write(*,*)'FINISHED: TLM_pbc_save_field'
 
  return
  
END SUBROUTINE TLM_pbc_save_field

!
! SUBROUTINE TLM_pbc_copy_saved_field
!
! Name TLM_pbc_copy_saved_field
!     
!
! Description
!     copy saved mesh data back into grid
!
! Comments:
!      
!
! History
!
!     started 29/4/14 CJS based on the original HIRF-SE code Fieldsolve
!     

SUBROUTINE TLM_pbc_copy_saved_field()

USE mesh
USE constants
USE TLM_general
USE TLM_periodic
USE cell_parameters
USE TLM_volume_materials
USE TLM_surface_materials

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer ix,iy,iz
  integer :: cx,cy,cz,face
  integer :: point,point2,mat,i

! START

!  write(*,*)'CALLED: TLM_pbc_copy_saved_field'
    
! volume stuff
    V(1:12,1:nx,1:ny,nz1:nz2)=V_pbc_save(1:12,1:nx,1:ny,nz1:nz2)

! boundary stuff
    V_pbc_x(1:ny,nz1:nz2,1:2,1:n_xbc_point)=V_pbc_x_save(1:ny,nz1:nz2,1:2,1:n_xbc_point)
    V_pbc_y(1:nx,nz1:nz2,1:2,1:n_ybc_point)=V_pbc_y_save(1:nx,nz1:nz2,1:2,1:n_ybc_point)
    
    xbc_write_point_ymin=xbc_write_point_save_ymin
    xbc_write_point_ymax=xbc_write_point_save_ymax
    
    ybc_write_point_xmin=ybc_write_point_save_xmin
    ybc_write_point_xmax=ybc_write_point_save_xmax


! need to copy surface material stuff here
	  
! need to copy volume material stuff here
        
!  write(*,*)'FINISHED: TLM_pbc_copy_saved_field'
 
 RETURN

END SUBROUTINE TLM_pbc_copy_saved_field
