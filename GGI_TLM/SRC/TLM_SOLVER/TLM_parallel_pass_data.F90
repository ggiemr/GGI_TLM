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
!SUBROUTINE TLM_parallel_pass_data_1()
!SUBROUTINE TLM_parallel_pass_data_2()
!
! Name TLM_parallel_pass_data_1
!    
!
! Description
!     Pass TLM link line voltage pulses from rank n to rank n+1 process
!     This is called after Scatter but before connect
!
! Comments:
!      
!
! History
!
!     started 08/07/09 CJS
!     Adapted for GGI_TLM 23/11/2012 CJS
!     Split parallel_pass_data into the more efficient two stage process 13/12/2012 CJS
!

SUBROUTINE TLM_parallel_pass_data_1()

USE mesh
USE cell_parameters
USE constants
USE TLM_general

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer ix,iy,iz
  integer pt,pts,npts,npts2
  INTEGER talk_to,nreal

! START

  CALL write_line('CALLED: parallel_pass_data_1',0,timestepping_output_to_screen_flag)

#if defined(MPI)  

! put the voltages to be sent into a 1D array  
  
  npts=nx*ny
  npts2=2*npts
  nreal=npts2
  
! voltage pulses to be sent to rank+1 processor  
  if (rank.ne.np-1) then
    do ix=1,nx
      do iy=1,ny
      
        pt=ix+(iy-1)*nx
    
        iz=nz2
        Vr_zmax(pt     )=V(Vx_zmax,ix,iy,iz)
        Vr_zmax(pt+npts)=V(Vy_zmax,ix,iy,iz)
   
      end do    ! next y cell
    end do      ! next x cell
  end if

! Even numbered processes talk to the rank+1 process, odd numbers to the rank-1 process in stage 1.    
  if (MOD(rank,2).EQ.0) then
    talk_to=rank+1
  else
    talk_to=rank-1
  end if
  
  if (talk_to.lt.np) THEN
  
    if (MOD(rank,2).EQ.0) then
! Even number process sends to rank+1 process

      CALL MPI_SEND(Vr_zmax,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)

    else  ! MOD(rank,2).NE.0
! Odd number process receives from rank-1 process

      CALL MPI_RECV(Vi_zmin,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
    
    end if ! MOD(rank.2)

  end if ! talk_to .lt. np    
   
! Even numbered processes talk to the rank-1 process, odd numbers to the rank+1 process in stage 2.    
  if (MOD(rank,2).EQ.0) then
    talk_to=rank-1
  else
    talk_to=rank+1
  end if
  
  if ( (talk_to.ge.0).AND.(talk_to.lt.np) ) THEN
  
    if (MOD(rank,2).EQ.0) then
! Even number process receives from rank-1 process

      CALL MPI_RECV(Vi_zmin,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
      
    else  ! MOD(rank,2).NE.0
! Odd number process sends to rank+1 process

      CALL MPI_SEND(Vr_zmax,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
   
    end if ! MOD(rank.2)

  end if ! talk_to .lt. np    
  
! voltage pulses received from rank-1 processor  
  if (rank.ne.0) then
    do ix=1,nx
      do iy=1,ny
      
        pt=ix+(iy-1)*nx
    
        iz=nz1-1
        V(Vx_zmax,ix,iy,iz)=Vi_zmin(pt     )
        V(Vy_zmax,ix,iy,iz)=Vi_zmin(pt+npts)
   
      end do    ! next y cell
    end do      ! next x cell
  end if

#endif

  CALL write_line('FINISHED: parallel_pass_data_1',0,timestepping_output_to_screen_flag)

  RETURN
  
END SUBROUTINE TLM_parallel_pass_data_1
!
! Name TLM_parallel_pass_data_2
!    
!
! Description
!     Pass TLM link line voltage pulses from rank n+1 to rank n process
!     This is called after Scatter but before connect
!
! Comments:
!      
!
! History
!
!     started 08/07/09 CJS
!     Adapted for GGI_TLM 23/11/2012 CJS
!     Split parallel_pass_data into the more efficient two stage process 13/12/2012 CJS
!

SUBROUTINE TLM_parallel_pass_data_2()

USE mesh
USE cell_parameters
USE constants
USE TLM_general

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer ix,iy,iz
  integer pt,pts,npts,npts2
  INTEGER talk_to,nreal

! START

  CALL write_line('CALLED: parallel_pass_data_2',0,timestepping_output_to_screen_flag)
  
#if defined(MPI)  
  
! put the voltages to be sent into a 1D array  
  
  npts=nx*ny
  npts2=2*npts
  nreal=npts2
    
! voltage pulses to be sent to rank-1 processor  
  if (rank.ne.0) then
    do ix=1,nx
      do iy=1,ny
      
        pt=ix+(iy-1)*nx
    
        iz=nz1-1
        Vr_zmin(pt     )=V(Vx_zmax,ix,iy,iz)
        Vr_zmin(pt+npts)=V(Vy_zmax,ix,iy,iz)
   
      end do    ! next y cell
    end do      ! next x cell
  end if

! Even numbered processes talk to the rank+1 process, odd numbers to the rank-1 process in stage 1.    
  if (MOD(rank,2).EQ.0) then
    talk_to=rank+1
  else
    talk_to=rank-1
  end if
  
  if (talk_to.lt.np) THEN
  
    if (MOD(rank,2).EQ.0) then
! Even number receives from the rank+1 process

      CALL MPI_RECV(Vi_zmax,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)      

    else  ! MOD(rank,2).NE.0
! Odd number process sends to the rank-1 process

      CALL MPI_SEND(Vr_zmin,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
    
    end if ! MOD(rank.2)

  end if ! talk_to .lt. np    
!   
! Even numbered processes talk to the rank-1 process, odd numbers to the rank+1 process in stage 2.    
  if (MOD(rank,2).EQ.0) then
    talk_to=rank-1
  else
    talk_to=rank+1
  end if
  
  if ( (talk_to.ge.0).AND.(talk_to.lt.np) ) THEN
  
    if (MOD(rank,2).EQ.0) then
! Even number process sends to the rank-1 process

      CALL MPI_SEND(Vr_zmin,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
      
    else  ! MOD(rank,2).NE.0
! Odd number process receives from the rank+1 process

      CALL MPI_RECV(Vi_zmax,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
   
    end if ! MOD(rank.2)

  end if ! talk_to .lt. np    
    
! voltage pulses received from rank+1 processor  
  if (rank.ne.np-1) then
    do ix=1,nx
      do iy=1,ny
      
        pt=ix+(iy-1)*nx
    
        iz=nz2
        V(Vx_zmax,ix,iy,iz)=Vi_zmax(pt     )
        V(Vy_zmax,ix,iy,iz)=Vi_zmax(pt+npts)
   
      end do    ! next y cell
    end do      ! next x cell
  end if

#endif

  CALL write_line('FINISHED: parallel_pass_data_2',0,timestepping_output_to_screen_flag)

  RETURN
  
END SUBROUTINE TLM_parallel_pass_data_2
