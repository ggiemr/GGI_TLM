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
! SUBROUTINE initialise_wrap_outer_boundary
! SUBROUTINE wrap_outer_boundary
!
! NAME
!     initialise_wrap_outer_boundary
!
! DESCRIPTION
!     Initialise wrapping boundary conditions on outer boundary
!     i.e. allocate memory for passing voltage pulse data for parallel runs
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 25/04/2012 CJS
!
!
SUBROUTINE initialise_wrap_outer_boundary

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

! START


! local_variables

  integer cx,cy,cz
  
  real*8 Vx,Vy,Vz

! START
  
  CALL write_line('CALLED: initialise_wrap_outer_boundary',0,timestepping_output_to_screen_flag)
  
  if (wrap_z) then  

! np.NE.1, we need to send data using mpi
    if (np.ne.1) then
        
      ALLOCATE( Vx_wrap_zmin_send(1:nx,1:ny) )
      ALLOCATE( Vy_wrap_zmin_send(1:nx,1:ny) )
      ALLOCATE( Vx_wrap_zmax_send(1:nx,1:ny) )
      ALLOCATE( Vy_wrap_zmax_send(1:nx,1:ny) ) 
        
      ALLOCATE( Vx_wrap_zmin_rcv(1:nx,1:ny) )
      ALLOCATE( Vy_wrap_zmin_rcv(1:nx,1:ny) )
      ALLOCATE( Vx_wrap_zmax_rcv(1:nx,1:ny) )
      ALLOCATE( Vy_wrap_zmax_rcv(1:nx,1:ny) ) 
    
    end if ! np.ne.1
    
  end if  
  

  CALL write_line('CALLED: initialise_wrap_outer_boundary',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_wrap_outer_boundary
!
!
! NAME
!     Wrap_outer_boundary
!
! DESCRIPTION
!     Wrapping boundary conditions on outer boundary
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 25/04/2012 CJS
!
!
SUBROUTINE wrap_outer_boundary

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

! START


! local_variables

  integer cx,cy,cz
  
  real*8 Vx,Vy,Vz
  
  integer n_data,talk_to
  
! START
  
  CALL write_line('CALLED: wrap_outer_boundary',0,timestepping_output_to_screen_flag)
    
  if (wrap_x) then  
  
! swap voltage pulses on xmin and xmax faces
    do cz=nz1,nz2
      do cy=1,ny
    
	Vy                 = V(Vy_xmin,1,cy,cz)
	V(Vy_xmin,1,cy,cz) = V(Vy_xmax,nx,cy,cz)
	V(Vy_xmax,nx,cy,cz)= Vy 
    
	Vz                 = V(Vz_xmin,1,cy,cz)
	V(Vz_xmin,1,cy,cz) = V(Vz_xmax,nx,cy,cz)
	V(Vz_xmax,nx,cy,cz)= Vz 
		
      end do    ! next y cell
    end do      ! next z cell
    
  end if ! wrap_x
  
  if (wrap_y) then  
  
! ymin and ymax faces
    do cz=nz1,nz2
      do cx=1,nx
       
	Vx                 = V(Vx_ymin,cx,1,cz)
	V(Vx_ymin,cx,1,cz) = V(Vx_ymax,cx,ny,cz)
	V(Vx_ymax,cx,ny,cz)= Vx 
    
	Vz                 = V(Vz_ymin,cx,1,cz)
	V(Vz_ymin,cx,1,cz) = V(Vz_ymax,cx,ny,cz)
	V(Vz_ymax,cx,ny,cz)= Vz 
    	  
      end do  ! next x cell
    end do      ! next z cell
    
  end if ! wrap_y
  
  if (wrap_z) then  

! np=1, simple wrap as for x and y directions
    if (np.eq.1) then
    
      do cy=1,ny
        do cx=1,nx
      
          Vx                 = V(Vx_zmin,cx,cy,1)
	  V(Vx_zmin,cx,cy,1) = V(Vx_zmax,cx,cy,nz)
	  V(Vx_zmax,cx,cy,nz)= Vx
      
          Vy                 = V(Vy_zmin,cx,cy,1)
	  V(Vy_zmin,cx,cy,1) = V(Vy_zmax,cx,cy,nz)
	  V(Vy_zmax,cx,cy,nz)= Vy
  	  
        end do  ! next x cell
      end do    ! next y cell

#if defined(MPI)  
    
    else
    
! parallel,have to package and send data here

      n_data=nx*ny
      
      if (rank.eq.0) then

! packege up boundary data
        do cy=1,ny
          do cx=1,nx
	
            Vx_wrap_zmin_send(cx,cy)=V(Vx_zmin,cx,cy,1)
            Vy_wrap_zmin_send(cx,cy)=V(Vy_zmin,cx,cy,1)
	    	  
          end do  ! next x cell
        end do    ! next y cell
      
! send data to rank np-1 process first          
        talk_to=np-1

        CALL MPI_SEND(Vx_wrap_zmin_send,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
        CALL MPI_SEND(Vy_wrap_zmin_send,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
      
! receive data data from rank np-1 process second          

        CALL MPI_RECV(Vx_wrap_zmin_rcv,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
        CALL MPI_RECV(Vy_wrap_zmin_rcv,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)

! put received zmax boundary data onto zmin face
        do cy=1,ny
          do cx=1,nx
	
            V(Vx_zmin,cx,cy,1)=Vx_wrap_zmin_rcv(cx,cy)
            V(Vy_zmin,cx,cy,1)=Vy_wrap_zmin_rcv(cx,cy)
	    	  
          end do  ! next x cell
        end do    ! next y cell

      else if (rank.eq.(np-1)) then

! packege up boundary data
      
        do cy=1,ny
          do cx=1,nx
	
            Vx_wrap_zmax_send(cx,cy)=V(Vx_zmax,cx,cy,nz)
            Vy_wrap_zmax_send(cx,cy)=V(Vy_zmax,cx,cy,nz)
	    	  
          end do  ! next x cell
        end do    ! next y cell

        talk_to=0
! receive data from rank np-1 process first          

        CALL MPI_RECV(Vx_wrap_zmax_rcv,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
        CALL MPI_RECV(Vy_wrap_zmax_rcv,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)

! send data to rank np-1 process second          

        CALL MPI_SEND(Vx_wrap_zmax_send,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
        CALL MPI_SEND(Vy_wrap_zmax_send,n_data,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)

! put received zmax boundary data onto zmin face
        do cy=1,ny
          do cx=1,nx
	
            V(Vx_zmax,cx,cy,nz)=Vx_wrap_zmax_rcv(cx,cy)
            V(Vy_zmax,cx,cy,nz)=Vy_wrap_zmax_rcv(cx,cy)
	    	  
          end do  ! next x cell
        end do    ! next y cell

      
      end if  ! rank .eq.np-1 

#endif
          
    end if ! np.ne.1
    
  end if  
  

  CALL write_line('CALLED: wrap_outer_boundary',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE wrap_outer_boundary
