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
! SUBROUTINE set_frequency_output_surfaces_in_mesh
! SUBROUTINE initialise_frequency_output_surface
! SUBROUTINE face_output_frequency_output_surfaces
!
! NAME
!     set_frequency_output_surfaces_in_mesh
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
SUBROUTINE set_frequency_output_surfaces_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer		:: output_surface
  integer		:: surface_number
  integer		:: number_of_faces
  integer		:: output_face
  integer		:: cx,cy,cz,face
  type(cell_point)	:: output_face1
  type(cell_point)	:: output_face2

! START
  
  CALL write_line('CALLED: set_frequency_output_surfaces_in_mesh',0,output_to_screen_flag)
    
  if (n_frequency_output_surfaces.gt.0) then

! FREQUENCY OUTPUT SURFACES

    if (rank.eq.0) then
      write(info_file_unit,*)'Number of Frequency Output Surfaces=',n_frequency_output_surfaces
    end if
    
    do output_surface=1,n_frequency_output_surfaces
    
      surface_number=frequency_output_surface(output_surface)%surface_number
      number_of_faces=problem_surfaces(surface_number)%number_of_faces
      frequency_output_surface(output_surface)%number_of_faces=number_of_faces
      
      if (rank.eq.0) then
        write(info_file_unit,*)'Frequency Output Surface number',output_surface,' Number of cell_faces=',number_of_faces
      end if
      
! allocate face list for the output surface      
      ALLOCATE( frequency_output_surface(output_surface)%face_list(1:number_of_faces) )
      
      do output_face=1,number_of_faces
      
! copy the faces from the geometry structure to the local
! frequency_output_surface structure
	
        cx=problem_surfaces(surface_number)%face_list(output_face)%cell%i
        cy=problem_surfaces(surface_number)%face_list(output_face)%cell%j
        cz=problem_surfaces(surface_number)%face_list(output_face)%cell%k
        face=problem_surfaces(surface_number)%face_list(output_face)%point
 
! check the side of the surface on which we want output and change if required,
	
        if      (face.eq.face_xmin) then      
	  if (.NOT.frequency_output_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cx=cx-1
	    face=face_xmax
	  end if
         else if (face.eq.face_xmax) then	 
	  if (.NOT.frequency_output_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cx=cx+1
	    face=face_xmin
	  end if
        else if (face.eq.face_ymin) then	 
	  if (.NOT.frequency_output_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cy=cy-1
	    face=face_ymax
	  end if
        else if (face.eq.face_ymax) then	 
	  if (.NOT.frequency_output_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cy=cy+1
	    face=face_ymin
	  end if
        else if (face.eq.face_zmin) then	 
	  if (.NOT.frequency_output_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cz=cz-1
	    face=face_zmax
	  end if
        else if (face.eq.face_zmax) then	 
	  if (.NOT.frequency_output_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cz=cz+1
	    face=face_zmin
	  end if
        end if
  	     
! Set the output face in the local_surface_output array
! We must set the output point number on the min face 

        output_face1%cell%i=cx
        output_face1%cell%j=cy
        output_face1%cell%k=cz
        output_face1%point=face
	
        CALL get_min_face(output_face1,output_face2)
	
        cx=output_face2%cell%i
        cy=output_face2%cell%j
        cz=output_face2%cell%k
	face=output_face2%point

        if (rank.eq.cell_face_rank(cz,face)) then
! output point belongs to this processor
	  
          local_surface_output(cx  ,cy  ,cz  ,face)=1 
  
        end if ! output point belongs to this processor

! preserve the face which includes the normal information in the frequency_output_surface%face_list	
        frequency_output_surface(output_surface)%face_list(output_face)=output_face1
                 
      end do !next cell face in this surface	  
  
    end do ! next frequency_output surface

  end if ! n_frequency_output_surface.GT.0 
  
  CALL write_line('FINISHED: set_frequency_output_surfaces_in_mesh',0,output_to_screen_flag)
  
END SUBROUTINE set_frequency_output_surfaces_in_mesh
!
! NAME
!     initialise_frequency_output_surfaces
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
!
!
SUBROUTINE initialise_frequency_output_surfaces

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

  integer	:: output_surface
  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: output_face
  integer	:: n_frames

  integer 	:: face_number
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4

! START
  
  if (n_frequency_output_surfaces.gt.0) then

! OUTPUT SURFACES
    do output_surface=1,n_frequency_output_surfaces
    
      number_of_faces=frequency_output_surface(output_surface)%number_of_faces
      
      write(*,*)'Frequency output surface number',output_surface,' number of cell_faces=',number_of_faces
      
      ALLOCATE( frequency_output_surface(output_surface)%face_output_field_number_list(1:number_of_faces) )
      ALLOCATE( frequency_output_surface(output_surface)%value(1:number_of_faces) )
      
      frequency_output_surface(output_surface)%value(1:number_of_faces)=(0d0,0d0)
      
      do output_face=1,number_of_faces
   
        cx  =frequency_output_surface(output_surface)%face_list(output_face)%cell%i
        cy  =frequency_output_surface(output_surface)%face_list(output_face)%cell%j
        cz  =frequency_output_surface(output_surface)%face_list(output_face)%cell%k
        face=frequency_output_surface(output_surface)%face_list(output_face)%point 
	
	if (rank.eq.cell_face_rank(cz,face)) then
   	     
          if      (face.eq.face_xmin) then      
            frequency_output_surface(output_surface)%face_output_field_number_list(output_face)=	&
	    			local_surface_output(cx  ,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_xmax) then	 
            frequency_output_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx+1,cy  ,cz  ,face_xmin)
          else if (face.eq.face_ymin) then	 
            frequency_output_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz  ,face_ymin)
          else if (face.eq.face_ymax) then	 
            frequency_output_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy+1,cz  ,face_ymin)
          else if (face.eq.face_zmin) then	 
            frequency_output_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz  ,face_zmin)
          else if (face.eq.face_zmax) then	 
            frequency_output_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz+1,face_zmin)
          end if
	
	end if ! cell belongs to this process
		  
      end do !next cell face in this surface	  
      
! set output timestep information
      CALL set_output_time_information(frequency_output_surface(output_surface)%specified_timestep_information,	&
                                       frequency_output_surface(output_surface)%first_timestep,   &
                                       frequency_output_surface(output_surface)%last_timestep,    &
                                       frequency_output_surface(output_surface)%timestep_interval,	    &
                                       frequency_output_surface(output_surface)%specified_time_information,	    &
                                       frequency_output_surface(output_surface)%first_time,	    &		  
                                       frequency_output_surface(output_surface)%last_time,	    &				  
                                       frequency_output_surface(output_surface)%time_interval,	&
                                       frequency_output_surface(output_surface)%number_of_output_timesteps )
  
    end do ! next output surface

    if (rank.eq.0) then
! rank 0 process only: write header for surface field outputs. Could maybe use write_surface_mesh_list_vtk...
  
      OPEN(unit=frequency_output_surface_unit,file=trim(problem_name)//frequency_output_surface_extn)

      write(frequency_output_surface_unit,'(A)')'# NUMBER OF FREQUENCY OUTPUT SURFACES:'
      write(frequency_output_surface_unit,'(I10)')n_frequency_output_surfaces

    end if ! rank 0 process

  end if ! n_frequency_output_surfaces.GT.0
  

  RETURN

END SUBROUTINE initialise_frequency_output_surfaces

!
! SUBROUTINE face_output_frequency_output_surfaces
!
! NAME
!     face_output_frequency_output_surfaces
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
!     16/12/2014  CJS  Fix bug in surface current density output - normal directions needed to be specified.  
!
!
SUBROUTINE face_output_frequency_output_surfaces

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer 	:: output_surface
  integer 	:: output_face
  integer 	:: output_field_number
  integer 	:: number_of_faces
  integer 	:: field_component
  integer 	:: face
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
  
  CALL write_line('CALLED: face_output_frequency_output_surfaces',0,timestepping_output_to_screen_flag)

  
  do output_surface=1,n_frequency_output_surfaces
  
    CALL get_output_flag(output_flag,	&
                         frequency_output_surface(output_surface)%first_timestep,	&
                         frequency_output_surface(output_surface)%last_timestep,	&
                         frequency_output_surface(output_surface)%timestep_interval )
			 
    if (output_flag) then
      
      number_of_faces=frequency_output_surface(output_surface)%number_of_faces
      frequency=frequency_output_surface(output_surface)%frequency
      ejwt=exp(-j*2d0*pi*frequency*time)
      
      do output_face=1,number_of_faces
         
        output_field_number=frequency_output_surface(output_surface)%face_output_field_number_list(output_face)
        field_component=frequency_output_surface(output_surface)%field_component  
        face=frequency_output_surface(output_surface)%face_list(output_face)%point
        cz  =frequency_output_surface(output_surface)%face_list(output_face)%cell%k
      
        if (rank.eq.cell_face_rank(cz,face)) then
! the output point is in this process so set the value   	     
	
          normx=0
	  normy=0
	  normz=0
	
          if	(face.eq.face_xmin) then
	    side=1
            normx=1
          else if (face.eq.face_xmax) then
	    side=2
            normx=-1
          else if (face.eq.face_ymin) then
	    side=1	
	    normy=1
          else if (face.eq.face_ymax) then
	    side=2	
	    normy=-1
          else if (face.eq.face_zmin) then
	    side=1	
	    normz=1
          else if (face.eq.face_zmax) then
	    side=2	
	    normz=-1
          end if
 
	  field(1:6)=face_output_field(output_field_number,side,1:6)
	
	  if (field_component.LE.6) then

            value=field(field_component)  

          else if (field_component.EQ.Jx) then
	
	    Jsx= ( normy*field(Hz)-normz*field(Hy) )   ! J=nxH
	    value=Jsx
	
          else if (field_component.EQ.Jy) then
	
	    Jsy= ( normz*field(Hx)-normx*field(Hz) )   ! J=nxH
	    value=Jsy
	
          else if (field_component.EQ.Jz) then
	
	    Jsz= ( normx*field(Hy)-normy*field(Hx) )   ! J=nxH
	    value=Jsz
		  
	  end if
	
	  frequency_output_surface(output_surface)%value(output_face)=	&
	      frequency_output_surface(output_surface)%value(output_face)+value*ejwt   !*dt

        end if ! output face belongs to this process
          
      end do !next cell face in this surface	  

    end if ! this timestep should contribute to the output
  
  end do ! next frequency_output_surface

  
  CALL write_line('FINISHED: face_output_frequency_output_surfaces',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_output_frequency_output_surfaces

! Name write_frequency_output_surfaces
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

SUBROUTINE write_frequency_output_surfaces

USE file_information
USE TLM_output
USE mesh
Use TLM_general
Use constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer :: output_surface
  
  integer :: output_face
  
  integer :: cx,cy,cz,face
  
  complex*16 :: value
 
  integer,allocatable	:: n_surfaces_rank(:)
  
  integer	:: talk_to
  integer	:: n_integer
  integer	:: loop
  integer 	:: surface_count
  integer 	:: surface_count_save
  
  type(xyz)	:: point1,point2,point3,point4
  integer 	:: point_number
  
  integer					:: local_n_faces
  type(frequency_output_surface_type)    	:: local_frequency_output_surface
  
  integer,allocatable	:: integer_array(:)
  complex*16,allocatable:: complex_array(:)

! START

  write(*,*)'CALLED: write_frequency_output_surfaces'

! write frequency_output_surface_data

  if (n_frequency_output_surfaces.ne.0) then

    do output_surface=1,n_frequency_output_surfaces
     
! Stage 1: send the number of surfaces in each process to rank 0
    
      if (rank.eq.0) then
      
        ALLOCATE ( n_surfaces_rank(0:np-1) )    
        n_surfaces_rank(0)=frequency_output_surface(output_surface)%number_of_faces
      
      end if
      
#if defined(MPI)     

      if (rank.ne.0) then 
    
        talk_to=0
        n_integer=1
        CALL MPI_SEND(frequency_output_surface(output_surface)%number_of_faces,n_integer,	&
	              MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

      else if (rank.eq.0) then 
      
        do talk_to=1,np-1
          n_integer=1
          CALL MPI_RECV(n_surfaces_rank(talk_to),n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
        end do
    
      end if
      
#endif

! Stage 2: Count the total number of frequency output surfaces
      if (rank.eq.0) then
        local_n_faces=0
        do loop=0,np-1
          local_n_faces=local_n_faces+n_surfaces_rank(loop)
        end do
	
! Stage 3: Allocate memory for the whole output dataset in the rank 0 process
        ALLOCATE( local_frequency_output_surface%face_list(1:local_n_faces) )
        ALLOCATE( local_frequency_output_surface%value(1:local_n_faces) )

      end if
      
! STAGE 4: communicate all of the data to the rank 0 process      
      if (rank.eq.0) then 
  
! set values already in rank 0 process    
        surface_count=0
        do loop=1,n_surfaces_rank(0)
          surface_count=surface_count+1
	  local_frequency_output_surface%face_list(surface_count)%cell%i=	&
	    frequency_output_surface(output_surface)%face_list(surface_count)%cell%i
	  local_frequency_output_surface%face_list(surface_count)%cell%j=	&
	    frequency_output_surface(output_surface)%face_list(surface_count)%cell%j
	  local_frequency_output_surface%face_list(surface_count)%cell%k=	&
	    frequency_output_surface(output_surface)%face_list(surface_count)%cell%k
	  local_frequency_output_surface%face_list(surface_count)%point=	&
	    frequency_output_surface(output_surface)%face_list(surface_count)%point
	  local_frequency_output_surface%value(surface_count)=	&
	    frequency_output_surface(output_surface)%value(surface_count)
        end do ! next rank 0 surface
	
      end if
      
#if defined(MPI)     

      if (rank.ne.0) then 

! send cx,cy,cz,face,value to rank 0 process    
        talk_to=0
        n_integer=frequency_output_surface(output_surface)%number_of_faces
        ALLOCATE( integer_array(1:n_integer) )
        ALLOCATE( complex_array(1:n_integer) )
! send cx      
        integer_array(1:n_integer)=frequency_output_surface(output_surface)%face_list(1:n_integer)%cell%i
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cy   
        integer_array(1:n_integer)=frequency_output_surface(output_surface)%face_list(1:n_integer)%cell%j
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cz    
        integer_array(1:n_integer)=frequency_output_surface(output_surface)%face_list(1:n_integer)%cell%k
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send face   
        integer_array(1:n_integer)=frequency_output_surface(output_surface)%face_list(1:n_integer)%point
        CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send value  
        complex_array(1:n_integer)=frequency_output_surface(output_surface)%value(1:n_integer)
        CALL MPI_SEND(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

        DEALLOCATE( integer_array )
        DEALLOCATE( complex_array )
      
      else if (rank.eq.0) then 
        
        do talk_to=1,np-1
      
          surface_count_save=surface_count
	  
          n_integer=n_surfaces_rank(talk_to)
          ALLOCATE( integer_array(1:n_integer) )
          ALLOCATE( complex_array(1:n_integer) )
! get cx
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          surface_count=surface_count_save
          do loop=1,n_integer
            surface_count=surface_count+1
	    local_frequency_output_surface%face_list(surface_count)%cell%i=integer_array(loop)
          end do
! get cy
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          surface_count=surface_count_save
          do loop=1,n_integer
            surface_count=surface_count+1
	    local_frequency_output_surface%face_list(surface_count)%cell%j=integer_array(loop)
          end do
! get cz
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          surface_count=surface_count_save
          do loop=1,n_integer
            surface_count=surface_count+1
	    local_frequency_output_surface%face_list(surface_count)%cell%k=integer_array(loop)
          end do
! get face
          CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          surface_count=surface_count_save
          do loop=1,n_integer
            surface_count=surface_count+1
	    local_frequency_output_surface%face_list(surface_count)%point=integer_array(loop)
          end do
! get value
          CALL MPI_RECV(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	
          surface_count=surface_count_save
          do loop=1,n_integer
            surface_count=surface_count+1
	    local_frequency_output_surface%value(surface_count)=complex_array(loop)
          end do

          DEALLOCATE( integer_array )
          DEALLOCATE( complex_array )
	
        end do ! next process to talk to

! check that we have the correct amount of data  

        if (surface_count.ne.local_n_faces) then
          write(*,*)'Error in write_frequency_output_surfaces'
          write(*,*)'Data counting error',surface_count
          write(*,*)'Surface count',surface_count
          write(*,*)'local_n_faces',local_n_faces
          STOP
        end if
      
      end if ! rank 0 process

#endif

! STAGE 5: The rank 0 process writes the data to file

      if (rank.eq.0) then

! Write header data
        write(frequency_output_surface_unit,'(A,I10)')'# START SURFACE FIELD FILE TEMPLATE, OUTPUT SURFACE NUMBER:',	&
	                                              output_surface
        write(frequency_output_surface_unit,'(A)')'# Number of points:'
        write(frequency_output_surface_unit,'(I10)')local_n_faces*4
        write(frequency_output_surface_unit,'(A)')'# Number of faces:'
        write(frequency_output_surface_unit,'(I10)')local_n_faces
        write(frequency_output_surface_unit,'(A)')'# Number of frames:'
        write(frequency_output_surface_unit,'(I10)')24

! Write face then point data

        do output_face=1,local_n_faces
    
     	  CALL get_cell_face_corner_coordinates(local_frequency_output_surface%face_list(output_face), &
     					        point1,point2,point3,point4)    

     	  write(frequency_output_surface_unit,8000)point1%x,point1%y,point1%z
     	  write(frequency_output_surface_unit,8000)point2%x,point2%y,point2%z
     	  write(frequency_output_surface_unit,8000)point3%x,point3%y,point3%z
     	  write(frequency_output_surface_unit,8000)point4%x,point4%y,point4%z

8000      format(3E14.5)
     			   
        end do! next face 
    
        point_number=0
        do output_face=1,local_n_faces
    
     	  write(frequency_output_surface_unit,8010)4,point_number,point_number+1,point_number+2,point_number+3
     	  point_number=point_number+4
  
8010 	  format(I3,4I8)
     
        end do ! next face

! Write values
    
        do output_face=1,local_n_faces
         
          cx  =local_frequency_output_surface%face_list(output_face)%cell%i
          cy  =local_frequency_output_surface%face_list(output_face)%cell%j
          cz  =local_frequency_output_surface%face_list(output_face)%cell%k
          face=local_frequency_output_surface%face_list(output_face)%point
	
	  value=local_frequency_output_surface%value(output_face)

  	  write(frequency_output_surface_unit,8020)cx,cy,cz,    &
	  	      dble(value ),     &
	              dimag(value ), &
          	      abs( value )
8020      format(3I8,3E16.8)
        
        end do !next cell face in this surface	  
  
        DEALLOCATE ( n_surfaces_rank )
        DEALLOCATE( local_frequency_output_surface%face_list )
        DEALLOCATE( local_frequency_output_surface%value )
   
      end if ! rank.eq.0
  
    end do ! next frequency_output_surface
  
  end if !(n_frequency_output_surfaces.ne.0) 

  write(*,*)'FINISHED: write_frequency_output_surfaces'
  return
  
END SUBROUTINE write_frequency_output_surfaces

