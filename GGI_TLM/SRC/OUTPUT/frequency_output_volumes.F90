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
! SUBROUTINE set_frequency_output_volumes_in_mesh
! SUBROUTINE initialise_frequency_output_volume
! SUBROUTINE cell_output_frequency_output_volumes
!
! NAME
!     set_frequency_output_volumes_in_mesh
!
! DESCRIPTION
!
!     
! COMMENTS
!
! HISTORY
!
!     started 5/12/2012 CJS
!
!
SUBROUTINE set_frequency_output_volumes_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer		:: output_volume
  integer		:: volume_number
  integer		:: number_of_cells
  integer		:: output_cell
  integer		:: cx,cy,cz
  type(ijk)		:: output_cell1
  type(ijk)		:: output_cell2
  
  integer		:: total_number_of_cells
  integer		:: cell_count
  integer		:: cxmin,cymin,czmin

! START
  
  CALL write_line('CALLED: set_frequency_output_volumes_in_mesh',0,output_to_screen_flag)
    
  if (n_frequency_output_volumes.gt.0) then

! FREQUENCY OUTPUT VOLUMES

    if (rank.eq.0) then
      write(info_file_unit,*)'Number of Frequency Output volumes=',n_frequency_output_volumes
    end if
    
    do output_volume=1,n_frequency_output_volumes
    
      volume_number=frequency_output_volume(output_volume)%volume_number
      
      if (frequency_output_volume_output_every.Eq.1) then
        number_of_cells=problem_volumes(volume_number)%number_of_cells
	total_number_of_cells=number_of_cells
      else
! we have to count the cells
        total_number_of_cells=problem_volumes(volume_number)%number_of_cells
	cxmin=nx
	cymin=ny
	czmin=nz
        do output_cell=1,total_number_of_cells
      
! copy the cells from the geometry structure to the local
! frequency_output_volume structure
	
          cx=problem_volumes(volume_number)%cell_list(output_cell)%cell%i
          cy=problem_volumes(volume_number)%cell_list(output_cell)%cell%j
          cz=problem_volumes(volume_number)%cell_list(output_cell)%cell%k
	  
	  cxmin=min(cxmin,cx)
	  cymin=min(cymin,cy)
	  czmin=min(czmin,cz)
	  
	end do

        cell_count=0
	
        do output_cell=1,total_number_of_cells
      
! copy the cells from the geometry structure to the local
! frequency_output_volume structure
	
          cx=problem_volumes(volume_number)%cell_list(output_cell)%cell%i
          cy=problem_volumes(volume_number)%cell_list(output_cell)%cell%j
          cz=problem_volumes(volume_number)%cell_list(output_cell)%cell%k
	  
	  if ( (MOD(cx-cxmin,frequency_output_volume_output_every).EQ.0).AND.   &
	       (MOD(cy-cymin,frequency_output_volume_output_every).EQ.0).AND.   &
	       (MOD(cz-czmin,frequency_output_volume_output_every).EQ.0) ) then
! add this cell to the output list

	     cell_count=cell_count+1
	     
	  end if
	end do
	 
	number_of_cells=cell_count
	
      end if
      
      write(*,*)'Total number of cells in volume=',total_number_of_cells
      write(*,*)'Number of output cells=',number_of_cells
      
      frequency_output_volume(output_volume)%number_of_cells=number_of_cells
      
      if (rank.eq.0) then
        write(info_file_unit,*)'Frequency Output volume number',output_volume,' Number of cells=',number_of_cells
      end if
      
! allocate cell list for the output volume      
      ALLOCATE( frequency_output_volume(output_volume)%cell_list(1:number_of_cells) )
      
      cell_count=0
      
      do output_cell=1,total_number_of_cells
      
! copy the cells from the geometry structure to the local
! frequency_output_volume structure
	
        cx=problem_volumes(volume_number)%cell_list(output_cell)%cell%i
        cy=problem_volumes(volume_number)%cell_list(output_cell)%cell%j
        cz=problem_volumes(volume_number)%cell_list(output_cell)%cell%k
	
	if ( (MOD(cx-cxmin,frequency_output_volume_output_every).EQ.0).AND.   &
	     (MOD(cy-cymin,frequency_output_volume_output_every).EQ.0).AND.   &
	     (MOD(cz-czmin,frequency_output_volume_output_every).EQ.0) ) then
! add this cell to the output list

	   cell_count=cell_count+1
 
          if (rank.eq.cell_rank(cz)) then
! output point belongs to this processor
	  
            local_cell_output(cx  ,cy  ,cz  )=1 
  
          end if ! output point belongs to this processor

! preserve the cell which includes the normal information in the frequency_output_volume%cell_list	
          frequency_output_volume(output_volume)%cell_list(cell_count)%i=cx
          frequency_output_volume(output_volume)%cell_list(cell_count)%j=cy
          frequency_output_volume(output_volume)%cell_list(cell_count)%k=cz
	     
	end if
                 
      end do !next cell cell in this volume	  
  
    end do ! next frequency_output volume
    
    write(*,*)'Number of cells set in this frequency output volume=',cell_count
    
    if (rank.eq.0) then
      write(info_file_unit,*)'__________________________________________________________________'
      write(info_file_unit,*)' '
    end if

  end if ! n_frequency_output_volume.GT.0 
  
  CALL write_line('FINISHED: set_frequency_output_volumes_in_mesh',0,output_to_screen_flag)
  
END SUBROUTINE set_frequency_output_volumes_in_mesh
!
! NAME
!     initialise_frequency_output_volumes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/12/2012 CJS
!     5/11/2015 CJS	Include file compression stuff
!
SUBROUTINE initialise_frequency_output_volumes

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
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4

! START
  
  if (n_frequency_output_volumes.gt.0) then

! OUTPUT VOLUMES
    do output_volume=1,n_frequency_output_volumes
    
      number_of_cells=frequency_output_volume(output_volume)%number_of_cells
      
      write(*,*)'Frequency output volume number',output_volume,' number of cell_cells=',number_of_cells
      
      ALLOCATE( frequency_output_volume(output_volume)%cell_output_field_number_list(1:number_of_cells) )
      ALLOCATE( frequency_output_volume(output_volume)%value(1:number_of_cells) )
      
      frequency_output_volume(output_volume)%value(1:number_of_cells)=(0d0,0d0)
      
      do output_cell=1,number_of_cells
   
        cx  =frequency_output_volume(output_volume)%cell_list(output_cell)%i
        cy  =frequency_output_volume(output_volume)%cell_list(output_cell)%j
        cz  =frequency_output_volume(output_volume)%cell_list(output_cell)%k
	
	if (rank.eq.cell_rank(cz)) then
   	     
          frequency_output_volume(output_volume)%cell_output_field_number_list(output_cell)=	&
	       local_cell_output(cx,cy,cz) 
	
	end if ! cell belongs to this process
		  
      end do !next cell cell in this volume	  
      
! set output timestep information
      CALL set_output_time_information(frequency_output_volume(output_volume)%specified_timestep_information,	&
                                       frequency_output_volume(output_volume)%first_timestep,   &
                                       frequency_output_volume(output_volume)%last_timestep,    &
                                       frequency_output_volume(output_volume)%timestep_interval,	    &
                                       frequency_output_volume(output_volume)%specified_time_information,	    &
                                       frequency_output_volume(output_volume)%first_time,	    &		  
                                       frequency_output_volume(output_volume)%last_time,	    &				  
                                       frequency_output_volume(output_volume)%time_interval,	&
                                       frequency_output_volume(output_volume)%number_of_output_timesteps )
  
    end do ! next output volume

    if (rank.eq.0) then
! rank 0 process only: write header for volume field outputs. Could maybe use write_volume_mesh_list_vtk...

      if (frequency_output_volume_format.EQ.frequency_output_volume_format_normal) then
  
        CALL open_output_file_write(frequency_output_volume_unit,	&
             trim(problem_name)//frequency_output_volume_extn,compress_output_files)

        write(frequency_output_volume_unit,'(A)')'# NUMBER OF FREQUENCY OUTPUT VOLUMES:'
        write(frequency_output_volume_unit,'(I10)')n_frequency_output_volumes
	
      end if ! normal output format
      
    end if ! rank 0 process

  end if ! n_frequency_output_volumes.GT.0
  

  RETURN

END SUBROUTINE initialise_frequency_output_volumes

!
! SUBROUTINE cell_output_frequency_output_volumes
!
! NAME
!     cell_output_frequency_output_volumes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE cell_output_frequency_output_volumes

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
  integer 	:: cz
  integer	:: side
  integer	:: normx,normy,normz
  real*8	:: value
  
  real*8	:: field(1:6)
  real*8	:: Jsx,Jsy,Jsz,Jmag,Emag,Hmag,P
  
  real*8	:: frequency
  complex*16	:: ejwt
  
  logical	:: output_flag

! START
  
  CALL write_line('CALLED: cell_output_frequency_output_volumes',0,timestepping_output_to_screen_flag)

  
  do output_volume=1,n_frequency_output_volumes
  
    CALL get_output_flag(output_flag,	&
                         frequency_output_volume(output_volume)%first_timestep,	&
                         frequency_output_volume(output_volume)%last_timestep,	&
                         frequency_output_volume(output_volume)%timestep_interval )
			 
    if (output_flag) then
      
      number_of_cells=frequency_output_volume(output_volume)%number_of_cells
      frequency=frequency_output_volume(output_volume)%frequency
      ejwt=exp(-j*2d0*pi*frequency*time)
    
      do output_cell=1,number_of_cells
         
        output_field_number=frequency_output_volume(output_volume)%cell_output_field_number_list(output_cell)
        field_component=frequency_output_volume(output_volume)%field_component  
        cz  =frequency_output_volume(output_volume)%cell_list(output_cell)%k
      
        if (rank.eq.cell_rank(cz)) then
 
	  value=cell_output_field(output_field_number,field_component)
	
	  frequency_output_volume(output_volume)%value(output_cell)=	&
	      frequency_output_volume(output_volume)%value(output_cell)+value*ejwt ! *dt : Note dt factor removed

        end if ! output cell belongs to this process
          
      end do !next cell cell in this volume	  

    end if ! this timestep should contribute to the output
  
  end do ! next frequency_output_volume

  
  CALL write_line('FINISHED: cell_output_frequency_output_volumes',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cell_output_frequency_output_volumes

! Name write_frequency_output_volumes
!     
!
! Description
!     
!
! Comments:
!      
!
! History
!
!     started 5/12/2012 CJS adapted from Fieldsolve
!

SUBROUTINE write_frequency_output_volumes

USE file_information
USE TLM_output
USE mesh
Use TLM_general
Use constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer :: output_volume
  
  integer :: output_cell
  
  integer :: cx,cy,cz,cell
  
  complex*16 :: value
 
  integer,allocatable	:: n_volumes_rank(:)
  
  integer	:: talk_to
  integer	:: n_integer
  integer	:: loop
  integer 	:: volume_count
  integer 	:: volume_count_save
  
  type(xyz)	:: point1,point2,point3,point4,point5,point6,point7,point8
  integer 	:: point_number
  
  integer					:: local_n_cells
  type(frequency_output_volume_type)    	:: local_frequency_output_volume
  
  integer,allocatable	:: integer_array(:)
  complex*16,allocatable:: complex_array(:)
  
  character(LEN=2) :: number_string,trim_number_string

! START

  write(*,*)'CALLED: write_frequency_output_volumes'

! write frequency_output_volume_data

  if (n_frequency_output_volumes.ne.0) then

    do output_volume=1,n_frequency_output_volumes
    
! Stage 1: send the number of volumes in each process to rank 0

      if (rank.eq.0) then
  
        ALLOCATE ( n_volumes_rank(0:np-1) )
    
        n_volumes_rank(0)=frequency_output_volume(output_volume)%number_of_cells
      
      end if

#if defined(MPI)

      if (rank.ne.0) then 
    
        talk_to=0
        n_integer=1
        CALL MPI_SEND(frequency_output_volume(output_volume)%number_of_cells,n_integer,	&
	              MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

      else if (rank.eq.0) then 
    
        do talk_to=1,np-1
          n_integer=1
          CALL MPI_RECV(n_volumes_rank(talk_to),n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
        end do
    
      end if
      
#endif

! Stage 2: Count the total number of frequency output volumes
      if (rank.eq.0) then
        local_n_cells=0
        do loop=0,np-1
          local_n_cells=local_n_cells+n_volumes_rank(loop)
        end do
	
! Stage 3: Allocate memory for the whole output dataset in the rank 0 process
        ALLOCATE( local_frequency_output_volume%cell_list(1:local_n_cells) )
        ALLOCATE( local_frequency_output_volume%value(1:local_n_cells) )

      end if

! STAGE 4: communicate all of the data to the rank 0 process
      if (rank.eq.0) then 
  
! set values already in rank 0 process    
        volume_count=0
        do loop=1,n_volumes_rank(0)
          volume_count=volume_count+1
	  local_frequency_output_volume%cell_list(volume_count)%i=	&
	    frequency_output_volume(output_volume)%cell_list(volume_count)%i
	  local_frequency_output_volume%cell_list(volume_count)%j=	&
	    frequency_output_volume(output_volume)%cell_list(volume_count)%j
	  local_frequency_output_volume%cell_list(volume_count)%k=	&
	    frequency_output_volume(output_volume)%cell_list(volume_count)%k
	  local_frequency_output_volume%value(volume_count)=	&
	    frequency_output_volume(output_volume)%value(volume_count)	    
        end do ! next rank 0 process
	
      end if

#if defined(MPI)      

      if (rank.ne.0) then 

! send cx,cy,cz,cell,value to rank 0 process    
        talk_to=0
        n_integer=frequency_output_volume(output_volume)%number_of_cells
        ALLOCATE( integer_array(1:n_integer) )
        ALLOCATE( complex_array(1:n_integer) )
! send cx      
        integer_array(1:n_integer)=frequency_output_volume(output_volume)%cell_list(1:n_integer)%i
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cy   
        integer_array(1:n_integer)=frequency_output_volume(output_volume)%cell_list(1:n_integer)%j
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cz    
        integer_array(1:n_integer)=frequency_output_volume(output_volume)%cell_list(1:n_integer)%k
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send value  
        complex_array(1:n_integer)=frequency_output_volume(output_volume)%value(1:n_integer)
        CALL MPI_SEND(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

        DEALLOCATE( integer_array )
        DEALLOCATE( complex_array )
      
      else if (rank.eq.0) then 
        
        do talk_to=1,np-1
      
          volume_count_save=volume_count
	  
          n_integer=n_volumes_rank(talk_to)
          ALLOCATE( integer_array(1:n_integer) )
          ALLOCATE( complex_array(1:n_integer) )
! get cx
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_frequency_output_volume%cell_list(volume_count)%i=integer_array(loop)
          end do
! get cy
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_frequency_output_volume%cell_list(volume_count)%j=integer_array(loop)
          end do
! get cz
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_frequency_output_volume%cell_list(volume_count)%k=integer_array(loop)
          end do
! get value
          CALL MPI_RECV(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          volume_count=volume_count_save
          do loop=1,n_integer
            volume_count=volume_count+1
	    local_frequency_output_volume%value(volume_count)=complex_array(loop)
          end do

          DEALLOCATE( integer_array )
          DEALLOCATE( complex_array )
	
        end do ! next process to talk to

! check that we have the correct amount of data  

        if (volume_count.ne.local_n_cells) then
          write(*,*)'Error in write_frequency_output_volumes'
          write(*,*)'Data counting error',volume_count
          write(*,*)'volume count',volume_count
          write(*,*)'local_n_cells',local_n_cells
          STOP
        end if
      
      end if ! rank 0 process

#endif

! STAGE 5: The rank 0 process writes the data to file

      if (rank.eq.0) then
      
       if (frequency_output_volume_format.EQ.frequency_output_volume_format_normal) then

! NORMAL FORMAT OUTPUT

! Write header data
        write(frequency_output_volume_unit,'(A,I10)')'# START volume FIELD FILE TEMPLATE, OUTPUT volume NUMBER:',	&
	                                              output_volume
        write(frequency_output_volume_unit,'(A)')'# Number of points:'
        write(frequency_output_volume_unit,'(I10)')local_n_cells*8
        write(frequency_output_volume_unit,'(A)')'# Number of cells:'
        write(frequency_output_volume_unit,'(I10)')local_n_cells*6
        write(frequency_output_volume_unit,'(A)')'# Number of frames:'
        write(frequency_output_volume_unit,'(I10)')24

! Write cell then point data

        do output_cell=1,local_n_cells
    
     	  CALL get_cell_corner_coordinates(local_frequency_output_volume%cell_list(output_cell), &
     					   point1,point2,point3,point4,point5,point6,point7,point8)    

     	  write(frequency_output_volume_unit,8000)point1%x,point1%y,point1%z
     	  write(frequency_output_volume_unit,8000)point2%x,point2%y,point2%z
     	  write(frequency_output_volume_unit,8000)point3%x,point3%y,point3%z
     	  write(frequency_output_volume_unit,8000)point4%x,point4%y,point4%z
     	  write(frequency_output_volume_unit,8000)point5%x,point5%y,point5%z
     	  write(frequency_output_volume_unit,8000)point6%x,point6%y,point6%z
     	  write(frequency_output_volume_unit,8000)point7%x,point7%y,point7%z
     	  write(frequency_output_volume_unit,8000)point8%x,point7%y,point8%z

8000      format(3E14.5)
     			   
        end do! next cell 
    
        point_number=0
        do output_cell=1,local_n_cells
    
          write(frequency_output_volume_unit,8010)4,point_number  ,point_number+1,point_number+2,point_number+3  ! xmin face
          write(frequency_output_volume_unit,8010)4,point_number+4,point_number+7,point_number+6,point_number+5  ! xmax face
          write(frequency_output_volume_unit,8010)4,point_number  ,point_number+4,point_number+5,point_number+1  ! ymin face
          write(frequency_output_volume_unit,8010)4,point_number+3,point_number+2,point_number+6,point_number+7  ! ymax face
          write(frequency_output_volume_unit,8010)4,point_number  ,point_number+3,point_number+7,point_number+4  ! zmin face
          write(frequency_output_volume_unit,8010)4,point_number+1,point_number+5,point_number+6,point_number+2  ! zmax face
          point_number=point_number+8
  
8010 	  format(I3,4I12)
     
        end do ! next cell

! Write values
    
        do output_cell=1,local_n_cells
         
          cx  =local_frequency_output_volume%cell_list(output_cell)%i
          cy  =local_frequency_output_volume%cell_list(output_cell)%j
          cz  =local_frequency_output_volume%cell_list(output_cell)%k
	
	  value=local_frequency_output_volume%value(output_cell)

  	  write(frequency_output_volume_unit,8020)cx,cy,cz,    &
	  	      dble(value ),     &
	              dimag(value ), &
          	      abs( value )
8020      format(3I8,3E16.8)
        
        end do !next cell cell in this volume	  
      
       else if (frequency_output_volume_format.EQ.frequency_output_volume_format_xyz_field) then

! XYZ FORMAT OUTPUT

        write(number_string,'(I2)')output_volume
	trim_number_string=ADJUSTL(number_string)
        OPEN(unit=frequency_output_volume_unit,file=trim(problem_name)//frequency_output_volume_extn//trim(trim_number_string))

        do output_cell=1,local_n_cells
	
          CALL get_cell_centre_coordinate(local_frequency_output_volume%cell_list(output_cell),point1)
	
	  value=local_frequency_output_volume%value(output_cell)

  	  write(frequency_output_volume_unit,8030)point1%x,point1%y,point1%z,    &
	  	      dble(value ),     &
	              dimag(value ), &
          	      abs( value )
8030      format(6E14.4)
     			   
        end do! next cell 
	
	CLOSE(unit=frequency_output_volume_unit)
  
       end if
  
       DEALLOCATE ( n_volumes_rank )
       DEALLOCATE( local_frequency_output_volume%cell_list )
       DEALLOCATE( local_frequency_output_volume%value )
   
      end if ! rank.eq.0
  
    end do ! next frequency_output_volume
  
  end if !(n_frequency_output_volumes.ne.0) 

  write(*,*)'FINISHED: write_frequency_output_volumes'
  return
  
END SUBROUTINE write_frequency_output_volumes

