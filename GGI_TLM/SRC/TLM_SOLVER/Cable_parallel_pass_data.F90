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
!SUBROUTINE Cable_parallel_pass_data_1
!SUBROUTINE Cable_parallel_pass_data_2
!
! Name Cable_parallel_pass_data_1
!     
!
! Description
!    
!   Pass wire cell data at the edge of each processors mesh
!   
!   
! Comments:
! 
!
! History
!
!     started 10/07/09 CJS
!     adapted for GGI_TLM 27/11/2012
!     Split into two processes for more efficient parallel implementation 13/12/2012
!


SUBROUTINE Cable_parallel_pass_data_1

USE TLM_general
USE geometry
USE mesh
USE cell_parameters
USE file_information
USE constants
USE cables

IMPLICIT NONE

! variables passed to subroutine

! local_variables
       
   INTEGER talk_to,nreal
   INTEGER count,segment,local_cell,face,lpt,nw,i
      
! function_types
       
! START       

  CALL write_line('CALLED: Cable_parallel_pass_data_1',0,timestepping_output_to_screen_flag)

#if defined(MPI)

! collect voltage pulses to be sent to rank+1 processor  (wire_Vr_zmax)
  if (rank.ne.np-1) then
  
    lpt=0
    do count=1,n_zmax_segments_send
      segment=zmax_segment_list_send(count)
      nw=bundle_segment_list(segment)%n_conductors
!      write(*,*)'rank',rank,' Sending segment ',segment
      do i=1,nw
        lpt=lpt+1
        wire_Vr_zmax(lpt)=bundle_segment_list(segment)%Vlink(i)
      end do 
    end do
    
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

       CALL MPI_SEND(wire_Vr_zmax,n_zmax_reals_send,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
!
    else  ! MOD(rank,2).NE.0
! Odd number process receives from rank-1 process
 
      CALL MPI_RECV(wire_Vi_zmin,n_zmin_reals_rcv,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
!      write(*,8000)'p1 rank',rank,' gets from rank',talk_to,' wire_Vi_zmin(1)',wire_Vi_zmin(1),n_zmin_reals_rcv
    
    end if ! MOD(rank.2)

  end if ! talk_to .lt. np    
  
! Even numbered processes talk to the rank-1 process, odd numbers to the rank+1 process in stage 2.    
  if (MOD(rank,2)==0) then
    talk_to=rank-1
  else
    talk_to=rank+1
  end if
  
  if ( (talk_to.ge.0).AND.(talk_to.lt.np) ) THEN
  
    if (MOD(rank,2).EQ.0) then
! Even number process receives from rank-1 process

      CALL MPI_RECV(wire_Vi_zmin,n_zmin_reals_rcv,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
      
    else  ! MOD(rank,2).NE.0
! Odd number process sends to rank+1 process

      CALL MPI_SEND(wire_Vr_zmax,n_zmax_reals_send,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
  
    end if ! MOD(rank.2)

  end if ! talk_to .lt. np    
    
! voltage pulses received from rank-1 processor  (wire_Vi_zmin)
  if (rank.ne.0) then
  
    lpt=0
    do count=1,n_zmin_segments_rcv
      segment=zmin_segment_list_rcv(count)
      nw=bundle_segment_list(segment)%n_conductors
!      write(*,*)'rank',rank,' receiving segment ',segment
      do i=1,nw
        lpt=lpt+1
        bundle_segment_list(segment)%Vlink(i)=wire_Vi_zmin(lpt)
      end do 
    end do

  end if

#endif
  
  CALL write_line('FINISHED: Cable_parallel_pass_data_1',0,timestepping_output_to_screen_flag)

  RETURN
  
END SUBROUTINE Cable_parallel_pass_data_1

!
! Name Cable_parallel_pass_data_2
!     
!
! Description
!    
!   Pass wire cell data at the edge of each processors mesh
!   
!   
! Comments:
! 
!
! History
!
!     started 10/07/09 CJS
!     adapted for GGI_TLM 27/11/2012
!     Split into two processes for more efficient parallel implementation 13/12/2012
!


SUBROUTINE Cable_parallel_pass_data_2

USE TLM_general
USE geometry
USE mesh
USE cell_parameters
USE file_information
USE constants
USE cables

IMPLICIT NONE

! variables passed to subroutine

! local_variables
       
   INTEGER talk_to,nreal
   INTEGER count,segment,local_cell,face,lpt,nw,i
      
! function_types
       
! START       

  CALL write_line('CALLED: Cable_parallel_pass_data_2',0,timestepping_output_to_screen_flag)

#if defined(MPI)
  
! voltage pulses to be sent to rank-1 processor  (wire_Vr_zmin)
  if (rank.ne.0) then
  
    lpt=0
    do count=1,n_zmin_segments_send
      segment=zmin_segment_list_send(count)
      nw=bundle_segment_list(segment)%n_conductors
!      write(*,*)'rank',rank,' Sending segment ',segment
      do i=1,nw
        lpt=lpt+1
        wire_Vr_zmin(lpt)=bundle_segment_list(segment)%Vlink(i)
      end do 
    end do

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

      CALL MPI_RECV(wire_Vi_zmax,n_zmax_reals_rcv,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)      

    else  ! MOD(rank,2).NE.0
! Odd number process sends to the rank-1 process

      CALL MPI_SEND(wire_Vr_zmin,n_zmin_reals_send,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
    
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

      CALL MPI_SEND(wire_Vr_zmin,n_zmin_reals_send,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
      
    else  ! MOD(rank,2).NE.0
! Odd number process receives from the rank+1 process

      CALL MPI_RECV(wire_Vi_zmax,n_zmax_reals_rcv,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
   
    end if ! MOD(rank.2)

  end if ! talk_to .lt. np    
    
! voltage pulses received from rank+1 processor  
  if (rank.ne.np-1) then
  
    lpt=0
    do count=1,n_zmax_segments_rcv
      segment=zmax_segment_list_rcv(count)
      nw=bundle_segment_list(segment)%n_conductors
!      write(*,*)'rank',rank,' receiving segment ',segment
      do i=1,nw
        lpt=lpt+1
        bundle_segment_list(segment)%Vlink(i)=wire_Vi_zmax(lpt)
      end do 
    end do
    
  end if

#endif
  
  CALL write_line('FINISHED: Cable_parallel_pass_data_2',0,timestepping_output_to_screen_flag)

  RETURN
  
END SUBROUTINE Cable_parallel_pass_data_2

