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
! SUBROUTINE run_TLM
!
! NAME
!     run_TLM
!
! DESCRIPTION
!     Main update loop for the TLM solution
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE run_TLM

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
  
  CALL write_line('CALLED: run_TLM',0,output_to_screen_flag)
  
  timestepping_output_to_screen_flag=.FALSE.
!  timestepping_output_to_screen_flag=.TRUE.

  if (rank.eq.0) then
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'#START OF MESH DESCRIPTION'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Cell size, dl=',dl,' m'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Nx=',nx 
    write(info_file_unit,*)'Ny=',ny 
    write(info_file_unit,*)'Nz=',nz 
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Total number of cells=',nx*ny*nz
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Outer boundary reflection coefficients'
    write(info_file_unit,'(A7,F8.4,A7,F8.4)')'R_xmin=',R_xmin,' R_xmax=',R_xmax
    write(info_file_unit,'(A7,F8.4,A7,F8.4)')'R_ymin=',R_ymin,' R_ymax=',R_ymax
    write(info_file_unit,'(A7,F8.4,A7,F8.4)')'R_zmin=',R_zmin,' R_zmax=',R_zmax
    write(info_file_unit,*)''
    write(info_file_unit,*)'#END OF MESH DESCRIPTION'
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'#START OF RUN INFORMATION'
    write(info_file_unit,*)''
    write(info_file_unit,*)'GGI_TLM version:',trim(GGI_TLM_version)
    write(info_file_unit,*)'GGI_TLM date   :',trim(GGI_TLM_date)
    write(info_file_unit,*)'GGI_TLM compilation date:',trim(GGI_TLM_compilation_date)
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Timestep=',dt,' seconds'
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Number of timesteps=',n_timesteps
    write(info_file_unit,*)''
    write(info_file_unit,*)'Number of processes= ',np
    write(info_file_unit,*)''
    write(info_file_unit,*)'bicubic_warp_flag= ',bicubic_warp_flag
    write(info_file_unit,*)'frequency_scale_flag= ',frequency_scale_flag
    write(info_file_unit,*)'frequency_scale= ',frequency_scale
    write(info_file_unit,*)''
    if (run_ngspice) then
      write(info_file_unit,*)'Include link to ngspice'
      write(info_file_unit,*)'ngspice_timestep_reduction_factor=',ngspice_timestep_factor
      write(info_file_unit,*)''
    end if
    write(info_file_unit,*)'GGI_TLM solution Started:'
    call write_date_and_time(info_file_unit)
    write(info_file_unit,*)'' 
    write(*,*)'Problem name:',trim(problem_name)
    write(*,*)''
    write(*,*)'TLM solution Started:'
    call write_date_and_time(0)
    write(*,*)'' 
    
  end if
  
  if ( run_ngspice )  then    
  
    if (np.gt.1) then 	 
      write(*,*)'We cannot run GGI_TLM with ngspice in parallel at the moment...'
      STOP 1
    end if
    
#if defined(INCLUDE_NGSPICE)    

    CALL initialise_ngspice()
      
    CALL ngspice_init_connect()
    
#else
    write(*,*)'ERROR: GGI_TLM has not been compiled with ngspice.'
    STOP 1
#endif
    
  end if

#if defined(MPI)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
#endif

  if (rank.eq.0) then
  
#if defined(SEQ)
!    call system("ps -o pid,%cpu,%mem,rss,vsz,comm -C GGI_TLM_SEQ > GGI_TLM_memory_usage.txt ")
    call system("ps -o pid,%cpu,%mem,rss,vsz,comm -C GGI_TLM_SEQ > GGI_TLM_memory_usage.txt ")
#elif defined(MPI)
!    call system("ps u -C GGI_TLM_MPI > GGI_TLM_memory_usage.txt ")
    call system("ps -o pid,%cpu,%mem,rss,vsz,comm -C GGI_TLM_MPI > GGI_TLM_memory_usage.txt ")
#endif
    
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

    close(scratch_file_unit)
    
    write(*,'(A)')""
    write(*,'(A)')"Note that RSS (resident set size) and VSZ (virtual memory size) are in Kilobytes"
    write(*,'(A)')""

    write(info_file_unit,'(A)')""
    write(info_file_unit,'(A)')""
    write(*,'(A)')""
    write(*,'(A)')'____________________________________________________'
    write(*,'(A)')""
    flush(info_file_unit)
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
    
    if ( (set_small_to_zero).AND.(mod(timestep,n_small_to_zero).EQ.0) ) then
      CALL TLM_set_small_to_zero()
    end if

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
    
#if defined(INCLUDE_NGSPICE)    
    
    if ( run_ngspice )  then    
    
! Transfer scattered voltage pulses from TLM nodes at time-dt/2d0, to Ngspice        
      CALL update_TLM_to_ngspice()
 
! Run ngspice up to a breakpoint set to simulate a period of a full TLM timestep i.e. at the next scatter time
      CALL update_ngspice(time+dt/2d0)
      
! Calculate the Ngspice voltage        
      CALL update_ngspice_to_TLM()
       
    end if

#endif

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
  
#if defined(SEQ)
!    call system("ps -o pid,%cpu,%mem,rss,vsz,comm -C GGI_TLM_SEQ > GGI_TLM_memory_usage.txt ")
    call system("ps -o pid,%cpu,%mem,rss,vsz,comm -C GGI_TLM_SEQ > GGI_TLM_memory_usage.txt ")
#elif defined(MPI)
!    call system("ps u -C GGI_TLM_MPI > GGI_TLM_memory_usage.txt ")
    call system("ps -o pid,%cpu,%mem,rss,vsz,comm -C GGI_TLM_MPI > GGI_TLM_memory_usage.txt ")
#endif
    
    write(info_file_unit,'(A)')""
    write(info_file_unit,'(A)')"Memory Usage on completion:"
    write(info_file_unit,'(A)')""

    write(*,'(A)')'____________________________________________________'
    write(*,'(A)')""
    write(*,'(A)')"Memory Usage:"
    write(*,'(A)')""
    
    open(unit=scratch_file_unit,file='GGI_TLM_memory_usage.txt')
15   CONTINUE
    read(scratch_file_unit,'(A256)',end=16)ipline
    write(info_file_unit,'(A)')trim(ipline)
    write(*,'(A)')trim(ipline)
    
    GOTO 15
    
16   CONTINUE
    write(*,'(A)')""
    write(*,'(A)')"Note that RSS (resident set size) and VSZ (virtual memory size) are in Kilobytes"
    write(*,'(A)')""

    write(info_file_unit,'(A)')""
    write(info_file_unit,'(A)')""
    write(*,'(A)')""
    write(*,'(A)')'____________________________________________________'
    write(*,'(A)')""
    flush(info_file_unit)
 
    write(info_file_unit,*)'#END OF RUN INFORMATION'
    write(info_file_unit,*)'' 
    
  end if
    
#if defined(INCLUDE_NGSPICE)    
  
  if ( run_ngspice )  then    
  
    CALL finish_ngspice()
          
  end if

#endif
  
  CALL write_line('FINISHED: run_TLM',0,output_to_screen_flag)

  RETURN

END SUBROUTINE run_TLM
