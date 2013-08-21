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
! SUBROUTINE set_output_volumes_in_mesh
! SUBROUTINE initialise_output_volumes
! SUBROUTINE cell_output_volumes
!
! NAME
!     set_output_volumes_in_mesh
!
! DESCRIPTION
!     loop over all the required outputs and flag all the output cells 
!     in the array local_cell_output(i,j,k) 
!     
! COMMENTS
!
! HISTORY
!
!     started 11/02/2013 CJS
!
!
SUBROUTINE set_output_volumes_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer	:: output_volume
  integer	:: volume_number
  integer	:: number_of_cells
  integer	:: cx,cy,cz
  integer	:: output_cell

! START
  
  CALL write_line('CALLED: set_output_volumes_in_mesh',0,output_to_screen_flag)
    
  if (n_output_volumes.gt.0) then

! OUTPUT volumeS
    do output_volume=1,n_output_volumes
    
      volume_number=output_volumes(output_volume)%volume_number
      number_of_cells=problem_volumes(volume_number)%number_of_cells
      output_volumes(output_volume)%number_of_cells=number_of_cells
      
! allocate cell list for the output volume      
      ALLOCATE( output_volumes(output_volume)%cell_list(1:number_of_cells) )

! allocate values      
      ALLOCATE( output_volumes(output_volume)%value(1:number_of_cells) )
      
      do output_cell=1,number_of_cells
      
! copy the cells from the geometry structure to the local
! output_volume structure
	
        cx=problem_volumes(volume_number)%cell_list(output_cell)%cell%i
        cy=problem_volumes(volume_number)%cell_list(output_cell)%cell%j
        cz=problem_volumes(volume_number)%cell_list(output_cell)%cell%k
 
        if (rank.eq.cell_rank(cz)) then
! output point belongs to this processor
	  
          local_cell_output(cx,cy,cz)=1 
  
        end if ! output point belongs to this processor
	
        output_volumes(output_volume)%cell_list(output_cell)=	&
	       problem_volumes(volume_number)%cell_list(output_cell)%cell
                 
      end do !next cell in this volume	  
  
    end do ! next output volume

  end if ! n_output_volumes.GT.0

  
END SUBROUTINE set_output_volumes_in_mesh
!
! NAME
!     initialise_output_volumes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/02/2013 CJS
!
!
SUBROUTINE initialise_output_volumes

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE Cables
USE file_information

IMPLICIT NONE

! local variables

  integer	:: output_volume
  integer	:: cx,cy,cz
  integer	:: number_of_cells
  integer	:: output_cell
  integer	:: n_frames

  integer 	:: cell_number
  integer	:: face
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4,point5,point6,point7,point8
 
  integer,allocatable	:: n_volumes_rank(:)
  
  integer	:: talk_to
  integer	:: n_integer
  integer	:: loop
  integer 	:: volume_count
  integer 	:: volume_count_save
 
  integer			:: local_n_cells
  type(output_volume_type)    	:: local_output_volume
  
  integer,allocatable		:: integer_array(:)

! START

  
  CALL write_line('CALLED: initialise_output_volumes',0,output_to_screen_flag)
  
  if (n_output_volumes.eq.0) RETURN

  if (rank.eq.0) then
! rank 0 process only: write header for volume field outputs. Could maybe use write_volume_mesh_list_vtk...
  
    OPEN(unit=volume_field_output_unit,file=trim(problem_name)//volume_field_output_extn)
     
    write(volume_field_output_unit,'(A)')'# NUMBER OF OUTPUT VOLUMES'
    write(volume_field_output_unit,*)n_output_volumes
     
  end if
      
! LOOP OVER OUTPUT volumeS SETTING MESH DATA
  do output_volume=1,n_output_volumes
  
    number_of_cells=output_volumes(output_volume)%number_of_cells
    
    ALLOCATE( output_volumes(output_volume)%cell_output_field_number_list(1:number_of_cells) )
    
    do output_cell=1,number_of_cells
  
      cx  =output_volumes(output_volume)%cell_list(output_cell)%i
      cy  =output_volumes(output_volume)%cell_list(output_cell)%j
      cz  =output_volumes(output_volume)%cell_list(output_cell)%k

      if (rank.eq.cell_rank(cz)) then
      
  	output_volumes(output_volume)%cell_output_field_number_list(output_cell)= &
        		    local_cell_output(cx,cy,cz) 
           
      end if ! cell belongs to this process
        	
    end do !next cell in this volume	
  
! set output timestep information
    CALL set_output_time_information(output_volumes(output_volume)%specified_timestep_information,  &
  				     output_volumes(output_volume)%first_timestep,   &
  				     output_volumes(output_volume)%last_timestep,    &
  				     output_volumes(output_volume)%timestep_interval, 	  &
  				     output_volumes(output_volume)%specified_time_information,	  &
  				     output_volumes(output_volume)%first_time,	  &		
  				     output_volumes(output_volume)%last_time, 	  &
  				     output_volumes(output_volume)%time_interval,   &
  				     output_volumes(output_volume)%number_of_output_timesteps )
  
  end do ! next output volume

! SET UP FILE HEADER INFORMATION

  do output_volume=1,n_output_volumes
  	   
! Stage 1: send the number of cells in each process to rank 0

    if (rank.eq.0) then
  
      ALLOCATE ( n_volumes_rank(0:np-1) )
    
      n_volumes_rank(0)=output_volumes(output_volume)%number_of_cells
    
    end if
    
#if defined(MPI)

    if (rank.ne.0) then 
    
      talk_to=0
      n_integer=1
      CALL MPI_SEND(output_volumes(output_volume)%number_of_cells,n_integer,        &
        	    MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror) 	  

    else if (rank.eq.0) then 
    
      do talk_to=1,np-1
    	n_integer=1
    	CALL MPI_RECV(n_volumes_rank(talk_to),n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)  
      end do
    
    end if

#endif

! Stage 2: Count the total number of output volumes
    if (rank.eq.0) then
      local_n_cells=0
      do loop=0,np-1
        local_n_cells=local_n_cells+n_volumes_rank(loop)
      end do
	
! Stage 3: Allocate memory for the whole output dataset in the rank 0 process
      ALLOCATE( local_output_volume%cell_list(1:local_n_cells) )

    end if

! STAGE 4: communicate all of the data to the rank 0 process
    if (rank.eq.0) then 
  
! set values already in rank 0 process    
      volume_count=0
      do loop=1,n_volumes_rank(0)
        volume_count=volume_count+1
	local_output_volume%cell_list(volume_count)%i=	&
	  output_volumes(output_volume)%cell_list(volume_count)%i
	local_output_volume%cell_list(volume_count)%j=	&
	  output_volumes(output_volume)%cell_list(volume_count)%j
        local_output_volume%cell_list(volume_count)%k=	&
	  output_volumes(output_volume)%cell_list(volume_count)%k
      end do ! next rank 0 process

    end if
    
#if defined(MPI)

    if (rank.ne.0) then 

! send cx,cy,cz,cell,value to rank 0 process    
      talk_to=0
      n_integer=output_volumes(output_volume)%number_of_cells
      ALLOCATE( integer_array(1:n_integer) )
! send cx      
      integer_array(1:n_integer)=output_volumes(output_volume)%cell_list(1:n_integer)%i
      CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cy   
      integer_array(1:n_integer)=output_volumes(output_volume)%cell_list(1:n_integer)%j
      CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cz    
      integer_array(1:n_integer)=output_volumes(output_volume)%cell_list(1:n_integer)%k
      CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

      DEALLOCATE( integer_array )
      
    else if (rank.eq.0) then 
     
      do talk_to=1,np-1
      
        volume_count_save=volume_count
	  
        n_integer=n_volumes_rank(talk_to)
        ALLOCATE( integer_array(1:n_integer) )
! get cx
        CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
        volume_count=volume_count_save
        do loop=1,n_integer
          volume_count=volume_count+1
	  local_output_volume%cell_list(volume_count)%i=integer_array(loop)
        end do
! get cy
        CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
        volume_count=volume_count_save
        do loop=1,n_integer
          volume_count=volume_count+1
	  local_output_volume%cell_list(volume_count)%j=integer_array(loop)
        end do
! get cz
        CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
        volume_count=volume_count_save
        do loop=1,n_integer
          volume_count=volume_count+1
	  local_output_volume%cell_list(volume_count)%k=integer_array(loop)
        end do

        DEALLOCATE( integer_array )
	
      end do ! next process to talk to

! check that we have the correct amount of data  

      if (volume_count.ne.local_n_cells) then
        write(*,*)'Error in write_output_volumes'
        write(*,*)'Data counting error',volume_count
        write(*,*)'volume count',volume_count
        write(*,*)'local_n_cells',local_n_cells
        STOP
      end if
      
    end if ! rank 0 process

#endif

! STAGE 5: The rank 0 process writes the data to file

    if (rank.eq.0) then
! rank 0 process only: write header for volume field outputs. Could maybe use write_volume_mesh_list_vtk...
    
      write(volume_field_output_unit,'(A,I5)')'# START volume FIELD FILE TEMPLATE, OUTPUT volume NUMBER:',output_volume

      n_frames=output_volumes(output_volume)%number_of_output_timesteps

      number_of_cells=local_n_cells

      write(volume_field_output_unit,'(A)')'# Number of points:'
      write(volume_field_output_unit,'(I10)')number_of_cells*24
      write(volume_field_output_unit,'(A)')'# Number of faces:'
      write(volume_field_output_unit,'(I10)')number_of_cells*6
      write(volume_field_output_unit,'(A)')'# Number of frames:'
      write(volume_field_output_unit,'(I10)')n_frames

      do cell_number=1,number_of_cells
    
     	CALL get_cell_corner_coordinates(local_output_volume%cell_list(cell_number), &
     					      point1,point2,point3,point4,point5,point6,point7,point8)   

     	write(volume_field_output_unit,8000)point1%x,point1%y,point1%z
     	write(volume_field_output_unit,8000)point2%x,point2%y,point2%z
     	write(volume_field_output_unit,8000)point3%x,point3%y,point3%z
     	write(volume_field_output_unit,8000)point4%x,point4%y,point4%z
	
     	write(volume_field_output_unit,8000)point5%x,point5%y,point5%z
     	write(volume_field_output_unit,8000)point8%x,point8%y,point8%z
     	write(volume_field_output_unit,8000)point7%x,point7%y,point7%z
     	write(volume_field_output_unit,8000)point6%x,point6%y,point6%z
 
     	write(volume_field_output_unit,8000)point1%x,point1%y,point1%z
     	write(volume_field_output_unit,8000)point5%x,point5%y,point5%z
     	write(volume_field_output_unit,8000)point6%x,point6%y,point6%z
     	write(volume_field_output_unit,8000)point2%x,point2%y,point2%z

     	write(volume_field_output_unit,8000)point4%x,point4%y,point4%z
     	write(volume_field_output_unit,8000)point3%x,point3%y,point3%z
     	write(volume_field_output_unit,8000)point7%x,point7%y,point7%z
     	write(volume_field_output_unit,8000)point8%x,point8%y,point8%z

     	write(volume_field_output_unit,8000)point1%x,point1%y,point1%z
     	write(volume_field_output_unit,8000)point4%x,point4%y,point4%z
     	write(volume_field_output_unit,8000)point8%x,point8%y,point8%z
     	write(volume_field_output_unit,8000)point5%x,point5%y,point5%z

     	write(volume_field_output_unit,8000)point2%x,point2%y,point2%z
     	write(volume_field_output_unit,8000)point6%x,point6%y,point6%z
     	write(volume_field_output_unit,8000)point7%x,point7%y,point7%z
     	write(volume_field_output_unit,8000)point3%x,point3%y,point3%z 
     			   
      end do! next cell 
     
8000  format(3E14.5)
    
      point_number=0
      do cell_number=1,number_of_cells*6
    
        write(volume_field_output_unit,8010)4,point_number  ,point_number+1,point_number+2,point_number+3  
        point_number=point_number+4
  
8010 	format(I3,4I8)
     
      end do ! next cell

      write(volume_field_output_unit,'(A,I5)')'# END VOLUME FIELD FILE TEMPLATE, OUTPUT VOLUME NUMBER:',output_volume

      DEALLOCATE ( n_volumes_rank )
      DEALLOCATE ( local_output_volume%cell_list )
      
      output_volumes(output_volume)%frame_number=0
      
    end if ! rank 0 process

  end do ! next output volume

  RETURN

END SUBROUTINE initialise_output_volumes
!
! SUBROUTINE cell_output_volumes
!
! NAME
!     cell_output_volumes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/03/2013 CJS
!
!
SUBROUTINE cell_output_volumes

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer 	:: output_volume
  integer 	:: output_cell
  integer 	:: output_field_number
  integer 	:: number_of_cells
  integer 	:: field_component
  integer 	:: cell
  integer 	:: cx,cy,cz
  integer	:: side
  real*8	:: value
    
  logical	:: output_flag

  integer 	:: i,p_rank
  integer 	:: nreal
 
  integer,allocatable	:: n_volumes_rank(:)
  
  integer	:: talk_to
  integer	:: n_integer
  integer	:: loop
  integer 	:: volume_count
  integer 	:: volume_count_save
 
  integer			:: local_n_cells
  type(output_volume_type)    	:: local_output_volume
  
  real*8,allocatable		:: integer_array(:)
  real*8,allocatable		:: real_array(:)
  
  real*8			:: field(1:6)
  real*8			:: Jsx,Jsy,Jsz,Jmag,Emag,Hmag
  real*8			:: normx,normy,normz

! START
  
  CALL write_line('CALLED: cell_output_volumes',0,timestepping_output_to_screen_flag)

  do output_volume=1,n_output_volumes
  
    CALL get_output_flag(output_flag,	&
                         output_volumes(output_volume)%first_timestep,	&
                         output_volumes(output_volume)%last_timestep,	&
                         output_volumes(output_volume)%timestep_interval )
			 
    if (output_flag) then

! loop over the output volume cells in this process and set the value
      do output_cell=1,output_volumes(output_volume)%number_of_cells
         
        output_field_number=output_volumes(output_volume)%cell_output_field_number_list(output_cell)
        field_component=output_volumes(output_volume)%field_component  
	
	field(1:6)=cell_output_field(output_field_number,1:6)

        output_volumes(output_volume)%value(output_cell)=field(field_component)  

      end do ! next cell
	   
! Stage 1: send the number of volumes in each process to rank 0

      if (rank.eq.0) then
  
        ALLOCATE ( n_volumes_rank(0:np-1) )
    
        n_volumes_rank(0)=output_volumes(output_volume)%number_of_cells
      
      end if

#if defined(MPI)

      if (rank.ne.0) then 
    
        talk_to=0
        n_integer=1
        CALL MPI_SEND(output_volumes(output_volume)%number_of_cells,n_integer,	&
	              MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

      else if (rank.eq.0) then 
    
        do talk_to=1,np-1
          n_integer=1
          CALL MPI_RECV(n_volumes_rank(talk_to),n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
        end do
    
      end if

#endif

! Stage 2: Count the total number of output volumes
      if (rank.eq.0) then
        local_n_cells=0
        do loop=0,np-1
          local_n_cells=local_n_cells+n_volumes_rank(loop)
        end do
	
! Stage 3: Allocate memory for the whole output dataset in the rank 0 process
        ALLOCATE( local_output_volume%cell_list(1:local_n_cells) )
        ALLOCATE( local_output_volume%value(1:local_n_cells) )

      end if

! STAGE 4: communicate all of the data to the rank 0 process
      if (rank.eq.0) then 
  
! set values already in rank 0 process    
        volume_count=0
        do loop=1,n_volumes_rank(0)
          volume_count=volume_count+1
	  local_output_volume%cell_list(volume_count)%i=	&
	    output_volumes(output_volume)%cell_list(volume_count)%i
	  local_output_volume%cell_list(volume_count)%j=	&
	    output_volumes(output_volume)%cell_list(volume_count)%j
	  local_output_volume%cell_list(volume_count)%k=	&
	    output_volumes(output_volume)%cell_list(volume_count)%k
	  local_output_volume%value(volume_count)=	&
	    output_volumes(output_volume)%value(volume_count)
        end do ! next rank 0 process
	
      end if
      
#if defined(MPI)	

      if (rank.ne.0) then 

! send cx,cy,cz,cell,value to rank 0 process    
        talk_to=0
        n_integer=output_volumes(output_volume)%number_of_cells
        ALLOCATE( integer_array(1:n_integer) )
        ALLOCATE( real_array(1:n_integer) )
! send cx      
        integer_array(1:n_integer)=output_volumes(output_volume)%cell_list(1:n_integer)%i
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cy   
        integer_array(1:n_integer)=output_volumes(output_volume)%cell_list(1:n_integer)%j
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cz    
        integer_array(1:n_integer)=output_volumes(output_volume)%cell_list(1:n_integer)%k
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send value  
        real_array(1:n_integer)=output_volumes(output_volume)%value(1:n_integer)
        CALL MPI_SEND(real_array,n_integer,MPI_DOUBLE,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

        DEALLOCATE( integer_array )
        DEALLOCATE( real_array )
      
      else if (rank.eq.0) then 
      
        do talk_to=1,np-1
      
          volume_count_save=volume_count
	  
          n_integer=n_volumes_rank(talk_to)
          ALLOCATE( integer_array(1:n_integer) )
          ALLOCATE( real_array(1:n_integer) )
! get cx
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_output_volume%cell_list(volume_count)%i=integer_array(loop)
          end do
! get cy
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_output_volume%cell_list(volume_count)%j=integer_array(loop)
          end do
! get cz
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_output_volume%cell_list(volume_count)%k=integer_array(loop)
          end do
! get value
          CALL MPI_RECV(real_array,n_integer,MPI_DOUBLE,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_output_volume%value(volume_count)=real_array(loop)
          end do

          DEALLOCATE( integer_array )
          DEALLOCATE( real_array )
	
        end do ! next process to talk to

! check that we have the correct amount of data  

        if (volume_count.ne.local_n_cells) then
          write(*,*)'Error in write_output_volumes'
          write(*,*)'Data counting error',volume_count
          write(*,*)'volume count',volume_count
          write(*,*)'local_n_cells',local_n_cells
          STOP
        end if
      
      end if ! rank 0 process

#endif

! STAGE 5: The rank 0 process writes the data to file

      if (rank.eq.0) then
      
        output_volumes(output_volume)%frame_number=output_volumes(output_volume)%frame_number+1 
	    
        write(volume_field_output_unit,'(A)')'# START OF VOLUME FIELD OUTPUT DATA'

        write(volume_field_output_unit,'(A)')'# OUTPUT VOLUME NUMBER'
        write(volume_field_output_unit,*)output_volume

        write(volume_field_output_unit,'(A)')'# FRAME NUMBER'
        write(volume_field_output_unit,*)output_volumes(output_volume)%frame_number
      
        write(volume_field_output_unit,'(A)')'# NUMBER OF cellS'
        write(volume_field_output_unit,*)local_n_cells*6
	
        write(volume_field_output_unit,'(A)')'# volume FIELD DATA'
    
        do output_cell=1,local_n_cells
         
          cx  =local_output_volume%cell_list(output_cell)%i
          cy  =local_output_volume%cell_list(output_cell)%j
          cz  =local_output_volume%cell_list(output_cell)%k
	
	  value=local_output_volume%value(output_cell)
	
	  do i=1,6
            write(volume_field_output_unit,'(E14.5)')value
          end do
	  
        end do !next cell cell in this volume	  
  
        DEALLOCATE ( n_volumes_rank )
        DEALLOCATE( local_output_volume%cell_list )
        DEALLOCATE( local_output_volume%value )
	
        write(volume_field_output_unit,'(A,I5)')'# END OF VOLUME FIELD OUTPUT DATA, OUTPUT volume NUMBER:',output_volume

      end if ! rank.eq.0
	  
    end if ! output_flag=TRUE
  
  end do ! next output_volume
  
  CALL write_line('FINISHED: cell_output_volumes',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cell_output_volumes

