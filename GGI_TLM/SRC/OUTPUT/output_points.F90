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
! SUBROUTINE set_output_points_in_mesh
! SUBROUTINE initialise_output_points
! SUBROUTINE face_output_point
! SUBROUTINE cell_output_point
!
! NAME
!     set_output_points_in_mesh
!
! DESCRIPTION
!     loop over all the required outputs and flag all the output cells/ faces 
!     in the arrays local_cell_output(i,j,k) or local_surface_output(i,j,k,face)
!     as required
!     
! COMMENTS
!
! HISTORY
!
!     started 14/08/2012 CJS
!     Parallel 23/11/2012 CJS
!     separate output types 5/12/2012 CJS
!
!
SUBROUTINE set_output_points_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer	:: output_point
  integer	:: cx,cy,cz,face
  type(cell_point)	:: output_face1

! START
  
  CALL write_line('CALLED: set_output_points_in_mesh',0,output_to_screen_flag)
  
  if (n_output_points.GT.0) then

    do output_point=1,n_output_points
  
      cx=output_points(output_point)%cell_point%cell%i
      cy=output_points(output_point)%cell_point%cell%j
      cz=output_points(output_point)%cell_point%cell%k
      face=output_points(output_point)%cell_point%point
  
      if (face.eq.centre) then
! cell centre output

        output_points(output_point)%rank =cell_rank(cz)

        if (rank.eq.output_points(output_point)%rank) then
! output point belongs to this processor
  
          local_cell_output(cx,cy,cz)=1
    
          write(*,*)'Setting cell centre output point',output_point
          write(*,*)'Cell Coordinates:',cx,cy,cz
  
        end if ! output point belongs to this processor
  
      else
! must be an output point on a face
! set the output point number on the min face 

        CALL get_min_face(output_points(output_point)%cell_point,output_face1)
	
        cx=output_face1%cell%i
        cy=output_face1%cell%j
        cz=output_face1%cell%k
	face=output_face1%point

        output_points(output_point)%rank =cell_face_rank(cz,face)

        if (rank.eq.output_points(output_point)%rank) then
! output point belongs to this processor
	  
          local_surface_output(cx,cy,cz,face)=1 
	  
          write(*,*)'Setting cell face output point',output_point
          write(*,*)'Coordinates:',cx,cy,cz,' face:',face
  
        end if ! output point belongs to this processor
  
      end if ! centre or face output
  
    end do ! next output point

  end if ! n_output_points>0
  
END SUBROUTINE set_output_points_in_mesh
!
! NAME
!     initialise_output_points
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
SUBROUTINE initialise_output_points

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

  integer	:: output_point
  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: output_face
  integer	:: n_frames

  integer 	:: face_number
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4

! START
  
  if (n_output_points.GT.0) then

    do output_point=1,n_output_points
  
      cx=output_points(output_point)%cell_point%cell%i
      cy=output_points(output_point)%cell_point%cell%j
      cz=output_points(output_point)%cell_point%cell%k
      face=output_points(output_point)%cell_point%point
      
      if (face.eq.centre) then
! cell centre output
        if (rank.eq.cell_rank(cz)) then
  
          output_points(output_point)%cell_output_field_number=local_cell_output(cx,cy,cz)
	  
        end if
	
      else 
      
        if ( rank.eq.cell_face_rank(cz,face) ) then
! must be an output point on a face
! set the output point number on the min face and give a -sign if it refers to the opposite side of the face
           
          if	(face.eq.face_xmin) then      
            output_points(output_point)%face_output_field_number=local_surface_output(cx  ,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_xmax) then       
            output_points(output_point)%face_output_field_number=local_surface_output(cx+1,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_ymin) then       
            output_points(output_point)%face_output_field_number=local_surface_output(cx  ,cy  ,cz  ,face_ymin) 
          else if (face.eq.face_ymax) then       
            output_points(output_point)%face_output_field_number=local_surface_output(cx  ,cy+1,cz  ,face_ymin) 
          else if (face.eq.face_zmin) then       
            output_points(output_point)%face_output_field_number=local_surface_output(cx  ,cy  ,cz  ,face_zmin) 
          else if (face.eq.face_zmax) then       
            output_points(output_point)%face_output_field_number=local_surface_output(cx  ,cy  ,cz+1,face_zmin) 
          end if
	    
        end if
  
      end if ! centre or face output
      
! set output timestep information
      CALL set_output_time_information(output_points(output_point)%specified_timestep_information,	&
                                       output_points(output_point)%first_timestep,   &
                                       output_points(output_point)%last_timestep,    &
                                       output_points(output_point)%timestep_interval,	    &
                                       output_points(output_point)%specified_time_information,	    &
                                       output_points(output_point)%first_time,	    &		  
                                       output_points(output_point)%last_time,	    &				  
                                       output_points(output_point)%time_interval,	&
                                       output_points(output_point)%number_of_output_timesteps )
      
    end do ! next output point
  
    if (rank.eq.0) then
    
      OPEN(unit=field_output_unit,file=trim(problem_name)//field_output_extn)
  
      CALL write_time_domain_header_data(field_output_unit,n_output_points,n_timesteps)
      
    end if
    
  end if ! n_output_points>0

  RETURN

END SUBROUTINE initialise_output_points
!
! SUBROUTINE face_output_points
!
! NAME
!     face_output_points
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
SUBROUTINE face_output_points

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer 	:: output_point
  integer 	:: output_face
  integer 	:: output_field_number
  integer 	:: number_of_faces
  integer 	:: field_component
  integer 	:: face
  integer 	:: cz
  integer	:: side
   real*8	:: value
  
  logical	:: output_flag

  integer 	:: i,p_rank
  integer 	:: talk_to,nreal

! START
  
  CALL write_line('CALLED: face_output_points',0,timestepping_output_to_screen_flag)

#if defined(MPI)

  if (rank.ne.0) then
! send output data to rank 0 processor
  
    talk_to=0
    nreal=1
  
    do output_point=1,n_output_points
  
      CALL get_output_flag(output_flag,	&
                           output_points(output_point)%first_timestep,	&
                           output_points(output_point)%last_timestep,	&
                           output_points(output_point)%timestep_interval )
			 
      if (output_flag) then

        p_rank=output_points(output_point)%rank
        if (p_rank.eq.rank) then
	  
          face=output_points(output_point)%cell_point%point
  
          if (face.ne.centre) then
! cell centre output: get the correct field value to output
            output_field_number=output_points(output_point)%face_output_field_number
            field_component=output_points(output_point)%field_component  
   	     
            if      (face.eq.face_xmin) then
	      side=1
            else if (face.eq.face_xmax) then
	      side=2
            else if (face.eq.face_ymin) then
	      side=1 	  
            else if (face.eq.face_ymax) then
	      side=2 	  
            else if (face.eq.face_zmin) then
	      side=1 	  
            else if (face.eq.face_zmax) then
	      side=2 	  
            end if
	
            output_points(output_point)%value=face_output_field(output_field_number,side,field_component)   
	
            CALL MPI_SEND(output_points(output_point)%value,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)	
!            write(*,8000)'rank',rank,' sends to rank ',talk_to,' value',output_points(output_point)%value,nreal
8000    format(A,I5,A,I5,A,E16.8,I10)

          end if ! face.NE.centre

        end if ! the output point is in this process
	
      end if ! output_flag.eq.TRUE.
      
    end do  ! next output point

  end if
      
#endif
      
  if (rank.eq.0) then
! this is the rank 0 process so collect output data from the other processes in order

    do output_point=1,n_output_points
    
      CALL get_output_flag(output_flag,	&
                           output_points(output_point)%first_timestep,	&
                           output_points(output_point)%last_timestep,	&
                           output_points(output_point)%timestep_interval )
			 
      if (output_flag) then
      
        p_rank=output_points(output_point)%rank
        talk_to=p_rank
        nreal=1
	  
        face=output_points(output_point)%cell_point%point
  
        if (face.ne.centre) then

#if defined(MPI)
      
          if (p_rank.ne.0) then
	
! get the value from another process
            CALL MPI_RECV(output_points(output_point)%value,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
!            write(*,8000)'rank',rank,' gets from rank ',talk_to,'  value',output_points(output_point)%value,nreal

          end if

#endif
	  
          if (p_rank.eq.0) then
	  
! set value in rank 0 process
! face output: get the correct field value to output

            output_field_number=output_points(output_point)%face_output_field_number
            field_component=output_points(output_point)%field_component  
   	     
            if      (face.eq.face_xmin) then
	      side=1
            else if (face.eq.face_xmax) then
	      side=2
            else if (face.eq.face_ymin) then
	      side=1 	  
            else if (face.eq.face_ymax) then
	      side=2 	  
            else if (face.eq.face_zmin) then
	      side=1 	  
            else if (face.eq.face_zmax) then
	      side=2 	  
            end if
	
            output_points(output_point)%value=face_output_field(output_field_number,side,field_component)   
	  
          end if ! p_rank.ne.0
	  
        end if ! face.NE.centre
      
      end if ! output flag = TRUE
      
    end do  ! next output point
    
! write complete output data to file    

    do output_point=1,n_output_points
  
      CALL get_output_flag(output_flag,	&
                           output_points(output_point)%first_timestep,	&
                           output_points(output_point)%last_timestep,	&
                           output_points(output_point)%timestep_interval )
			 
      if (output_flag) then
  
        face=output_points(output_point)%cell_point%point
  
        if (face.ne.centre) then
! cell centre output
  
          value=output_points(output_point)%value
          if ( abs(value).lt.1D-30 )value=0d0

          write(field_output_unit,time_domain_output_format)time,output_point,value
    
        end if
      
      end if ! output_flag
  
    end do ! next output 
      
  end if ! (rank.eq.0)
  
  CALL write_line('FINISHED: face_output_points',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_output_points
!
! SUBROUTINE cell_output_point
!
! NAME
!     cell_output_point
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
!     parallel 23/08/2012 CJS
!
!
SUBROUTINE cell_output_point

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats

IMPLICIT NONE

! local variables

  integer 	:: output_point
  integer 	:: output_field_number
  integer 	:: field_component
  integer 	:: face
  integer 	:: cz
  real*8	:: value
  
  logical	:: output_flag

  integer 	:: i,p_rank
  integer 	:: talk_to,nreal

! START
  
  CALL write_line('CALLED: cell_output_point',0,timestepping_output_to_screen_flag)

#if defined(MPI)

  if (rank.ne.0) then
! send output data to rank 0 processor
  
    talk_to=0
    nreal=1
  
    do output_point=1,n_output_points
  
      CALL get_output_flag(output_flag,	&
                           output_points(output_point)%first_timestep,	&
                           output_points(output_point)%last_timestep,	&
                           output_points(output_point)%timestep_interval )
			 
      if (output_flag) then

        p_rank=output_points(output_point)%rank
        if (p_rank.eq.rank) then
	  
          face=output_points(output_point)%cell_point%point
  
          if (face.eq.centre) then
! cell centre output: get the correct field value to output
            output_field_number=output_points(output_point)%cell_output_field_number
            field_component=output_points(output_point)%field_component  
            output_points(output_point)%value=cell_output_field(output_field_number,field_component)      
	
            CALL MPI_SEND(output_points(output_point)%value,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)	
    !        write(*,8000)'rank',rank,' sends to rank ',talk_to,' value',output_points(output_point)%value,nreal
          end if ! face=centre

8000    format(A,I5,A,I5,A,E16.8,I10)

        end if ! the output point is in this process
	
      end if ! output_flag.eq.TRUE.
      
    end do  ! next output point

  end if ! rank.ne.0

#endif  
      
  if (rank.eq.0) then
! this is the rank 0 process so collect output data from the other processes in order

    do output_point=1,n_output_points
    
      CALL get_output_flag(output_flag,	&
                           output_points(output_point)%first_timestep,	&
                           output_points(output_point)%last_timestep,	&
                           output_points(output_point)%timestep_interval )
			 
      if (output_flag) then
      
        p_rank=output_points(output_point)%rank
        talk_to=p_rank
        nreal=1
	  
        face=output_points(output_point)%cell_point%point
  
        if (face.eq.centre) then

#if defined(MPI) 
     
          if (p_rank.ne.0) then
	
! get the value from another process
            CALL MPI_RECV(output_points(output_point)%value,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
!            write(*,8000)'rank',rank,' gets from rank ',talk_to,'  value',output_points(output_point)%value,nreal

          end if
	  
#endif
          if(p_rank.eq.0) then
 
            output_field_number=output_points(output_point)%cell_output_field_number
            field_component=output_points(output_point)%field_component  
            output_points(output_point)%value=cell_output_field(output_field_number,field_component)      
	   
          end if ! p_rank.ne.0
	
        end if ! face=centre
	
      end if ! output_flag=TRUE
      
    end do  ! next output point
    
! write complete output data to file    

    do output_point=1,n_output_points
  
      CALL get_output_flag(output_flag,	&
                           output_points(output_point)%first_timestep,	&
                           output_points(output_point)%last_timestep,	&
                           output_points(output_point)%timestep_interval )
			 
      if (output_flag) then
  
        face=output_points(output_point)%cell_point%point
  
        if (face.eq.centre) then
! cell centre output
  
          value=output_points(output_point)%value
          if ( abs(value).lt.1D-30 )value=0d0
          write(field_output_unit,time_domain_output_format)time,output_point,value
    
        end if
      
      end if ! output_flag
  
    end do ! next output 
      
  end if ! (rank.eq.0)
  
  CALL write_line('FINISHED: cell_output_point',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cell_output_point

