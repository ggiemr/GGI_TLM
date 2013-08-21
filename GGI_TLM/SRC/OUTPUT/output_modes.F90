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
! SUBROUTINE set_mode_outputs_in_mesh
!
! NAME
!     set_mode_outputs_in_mesh
!
! DESCRIPTION
!     
!     
! COMMENTS
!  
!
! HISTORY
!
!     started 11/1/2013 CJS
!
!
SUBROUTINE set_mode_outputs_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer	:: output_mode
  integer	:: surface_number
  
  integer	:: surface_cxmin,surface_cxmax,surface_nx
  integer	:: surface_cymin,surface_cymax,surface_ny
  integer	:: surface_czmin,surface_czmax,surface_nz
  
  integer	:: n_surface_cells
  
  integer	:: mode_cxmin,mode_cxmax,mode_nx
  integer	:: mode_cymin,mode_cymax,mode_ny
  integer	:: mode_czmin,mode_czmax,mode_nz
  
  integer	:: n_mode_cells
  
  integer	:: offset_x,offset_y,offset_z
  
  integer	:: n_field_samples,sample
  
  integer	:: n_mode_samples
  
  integer 	:: face
  
  integer 	:: i
  
  integer	:: cx,cy,cz
  real*8	:: Field
  
  integer		:: coordinates(1:3)
  integer		:: n_data_cols
  real*8,allocatable	:: data_cols(:)
  integer		:: xcol,ycol,zcol,mode_col
  
  integer	:: output_face
  type(cell_point)	:: output_face1
  type(cell_point)	:: output_face2
  
! START
  
  CALL write_line('CALLED: set_mode_outputs_in_mesh',0,output_to_screen_flag)
  
  if (n_output_modes.GT.0) then
  
    do output_mode=1,n_output_modes
    
! get the coordinate range of the specified surface
      surface_number=output_mode_list(output_mode)%surface_number
      
      surface_cxmin=problem_surfaces(surface_number)%mesh_cell_xmin
      surface_cxmax=problem_surfaces(surface_number)%mesh_cell_xmax
      surface_nx=surface_cxmax-surface_cxmin+1
      
      surface_cymin=problem_surfaces(surface_number)%mesh_cell_ymin
      surface_cymax=problem_surfaces(surface_number)%mesh_cell_ymax
      surface_ny=surface_cymax-surface_cymin+1
      
      surface_czmin=problem_surfaces(surface_number)%mesh_cell_zmin
      surface_czmax=problem_surfaces(surface_number)%mesh_cell_zmax
      surface_nz=surface_czmax-surface_czmin+1

! open the mode field file      
      
      open(UNIT=mode_file_unit,					      &
      	   FILE=output_mode_list(output_mode)%mode_file_name,     &
      	   STATUS='old')
	   
! work out the number of field samples in the mode field file and the cell range in x, y and z
      n_field_samples=0
      
      mode_cxmin=nx
      mode_cxmax=0
      mode_cymin=ny
      mode_cymax=0
      mode_czmin=nz
      mode_czmax=0
      
      xcol=output_mode_list(output_mode)%xcol
      ycol=output_mode_list(output_mode)%ycol
      zcol=output_mode_list(output_mode)%zcol
      mode_col=output_mode_list(output_mode)%mode_col-3
      
      n_data_cols=output_mode_list(output_mode)%mode_col-3
      ALLOCATE( data_cols(1:n_data_cols) )

10      CONTINUE

        read(mode_file_unit,*,end=100)coordinates(1:3),data_cols(1:n_data_cols)
	
	cx=coordinates(xcol)
	cy=coordinates(ycol)
	cz=coordinates(zcol)
	Field=data_cols(mode_col)
	
	mode_cxmin=min(mode_cxmin,cx)
	mode_cxmax=max(mode_cxmax,cx)
	mode_cymin=min(mode_cymin,cy)
	mode_cymax=max(mode_cymax,cy)
	mode_czmin=min(mode_czmin,cz)
	mode_czmax=max(mode_czmax,cz)
	
	n_field_samples=n_field_samples+1
	
	GOTO 10

100   CONTINUE	

      mode_nx=mode_cxmax-mode_cxmin+1
      mode_ny=mode_cymax-mode_cymin+1
      mode_nz=mode_czmax-mode_czmin+1
            
! check that mode surface is a plane

      if ( (mode_nx.ne.1).AND.	&
           (mode_ny.ne.1).AND.	&
           (mode_nz.ne.1) ) then
	
	write(*,*)'Error in set_mode_outputs_in_mesh:'
	write(*,*)'Mode surface is not a plane'
	write(*,*)' mode_nx=',mode_nx
	write(*,*)' mode_ny=',mode_ny
	write(*,*)' mode_nz=',mode_nz
	STOP
	   
      end if
      
! check that output surface is a plane

      if ( (surface_nx.ne.1).AND.	&
           (surface_ny.ne.1).AND.	&
           (surface_nz.ne.1) ) then
	
	write(*,*)'Error in set_mode_outputs_in_mesh:'
	write(*,*)'mode output surface is not a plane'
	write(*,*)' surface_nx=',surface_nx
	write(*,*)' surface_ny=',surface_ny
	write(*,*)' surface_nz=',surface_nz
	STOP
	   
      end if
      
! check that the cell ranges for surface and mode are the same

      if ( (mode_nx.ne.surface_nx).OR.	&
           (mode_ny.ne.surface_ny).OR.	&
           (mode_nz.ne.surface_nz) ) then
	
	write(*,*)'Error in set_mode_outputs_in_mesh:'
	write(*,*)'Discrepancy between mode surface dimension and surface dimension'
	write(*,*)'surface_nx=',surface_nx,' mode_nx=',mode_nx
	write(*,*)'surface_ny=',surface_ny,' mode_ny=',mode_ny
	write(*,*)'surface_nz=',surface_nz,' mode_nz=',mode_nz
	STOP
	   
      end if

! calculate offset between surface and mode cells           
      offset_x=surface_cxmin-mode_cxmin
      offset_y=surface_cymin-mode_cymin
      offset_z=surface_czmin-mode_czmin
      
! calculate the number of mode output cells
      n_surface_cells=surface_nx*surface_ny*surface_nz
      n_mode_cells=mode_nx*mode_ny*mode_nz

! check consistency      
      if (n_surface_cells.ne.n_mode_cells) then
	
	write(*,*)'Error in set_mode_outputs_in_mesh:'
	write(*,*)'Discrepancy between number of mode surface cells and output surface cells'
	write(*,*)'n_surface_cells=',n_surface_cells,' n_mode_cells=',n_mode_cells
	STOP
	   
      end if

! Removed check - fails in parallel as problem_surfaces(surface_number)%number_of_faces
! is the number of faces belonging to this process only.
!
!! check the number of surface cells
!      if (n_mode_cells.ne.problem_surfaces(surface_number)%number_of_faces) then
!	
!	write(*,*)'Error in set_mode_outputs_in_mesh:'
!	write(*,*)'Discrepancy in number of output surface cells'
!	write(*,*)'n_mode_cells=',n_mode_cells
!	write(*,*)'problem_surfaces(surface_number)%number_of_faces=',	&
!	           problem_surfaces(surface_number)%number_of_faces
!	   
!      end if

! get the face number and check that this is consistent with the surface normal

      if (problem_surfaces(surface_number)%number_of_faces.ne.0) then

        face=problem_surfaces(surface_number)%face_list(1)%point
	
        do i=2,problem_surfaces(surface_number)%number_of_faces
          if (problem_surfaces(surface_number)%face_list(i)%point.ne.face) then
	
	    write(*,*)'Error in set_mode_outputs_in_mesh:'
	    write(*,*)'Not all the surface cell faces are the same'
	    write(*,*)'face(1)=',face
	    write(*,*)'i=',i
	    write(*,*)'face(i)=',problem_surfaces(surface_number)%face_list(i)%point
	    STOP
	  
	  end if
        end do ! next surface cell face

        if ( (face.eq.face_xmin).OR.(face.eq.face_xmax) ) then
          if (mode_nx.ne.1) then
	    write(*,*)'Error in set_mode_outputs_in_mesh:'
	    write(*,*)'face normal=x and mode_nx.ne.1'
	    STOP
	  end if
        else if ( (face.eq.face_ymin).OR.(face.eq.face_ymax) ) then
          if (mode_ny.ne.1) then
	    write(*,*)'Error in set_mode_outputs_in_mesh:'
	    write(*,*)'face normal=y and mode_ny.ne.1'
	    STOP
	  end if     
        else if ( (face.eq.face_zmin).OR.(face.eq.face_zmax) ) then
          if (mode_nz.ne.1) then
	    write(*,*)'Error in set_mode_outputs_in_mesh:'
	    write(*,*)'face normal=z and mode_nz.ne.1'
	    STOP
	  end if      
        end if

! MOVED DOWN      end if ! n_faces.ne.0

! re-read the mode file and work out how many samples are in this process
     
        rewind(mode_file_unit)
      
        n_mode_samples=0
      
        do sample=1,n_field_samples

          read(mode_file_unit,*,end=100)coordinates(1:3),data_cols(1:n_data_cols)
	
	  cx=coordinates(xcol)+offset_x
	  cy=coordinates(ycol)+offset_y
	  cz=coordinates(zcol)+offset_z
	  Field=data_cols(mode_col)
	
	  if (cell_face_rank(cz,face).eq.rank) then
	    n_mode_samples=n_mode_samples+1
	  end if
	
        end do

! allocate memory for the mode_output      

        output_mode_list(output_mode)%n_mode_samples=n_mode_samples
  
        ALLOCATE( output_mode_list(output_mode)%face_list(1:n_mode_samples) )
        ALLOCATE( output_mode_list(output_mode)%mode_field(1:n_mode_samples) )
        ALLOCATE( output_mode_list(output_mode)%face_output_field_number_list(1:n_mode_samples) )
      
! re-read the mode file data into the mode_output structure
     
        rewind(mode_file_unit)     
      
        n_mode_samples=0
      
        do sample=1,n_field_samples

          read(mode_file_unit,*,end=100)coordinates(1:3),data_cols(1:n_data_cols)
	
	  cx=coordinates(xcol)+offset_x
	  cy=coordinates(ycol)+offset_y
	  cz=coordinates(zcol)+offset_z
	  Field=data_cols(mode_col)
      	
	  if (cell_face_rank(cz,face).eq.rank) then
	    n_mode_samples=n_mode_samples+1
	    output_mode_list(output_mode)%face_list(n_mode_samples)%cell%i=cx
	    output_mode_list(output_mode)%face_list(n_mode_samples)%cell%j=cy
	    output_mode_list(output_mode)%face_list(n_mode_samples)%cell%k=cz
	    output_mode_list(output_mode)%face_list(n_mode_samples)%point=face
	    output_mode_list(output_mode)%mode_field(n_mode_samples)=Field
	  end if
	
        end do ! next sample
      
      else
! no face samples in this process

        output_mode_list(output_mode)%n_mode_samples= 0  
	   
      end if ! n_faces.ne.0     
      
      DEALLOCATE( data_cols )
    
      close(unit=mode_file_unit)

! set the mode output faces in the mesh     
 
      do output_face=1,output_mode_list(output_mode)%n_mode_samples
      
! copy the faces from the geometry structure to the local
! output_surface structure
	
        cx=output_mode_list(output_mode)%face_list(output_face)%cell%i
        cy=output_mode_list(output_mode)%face_list(output_face)%cell%j
        cz=output_mode_list(output_mode)%face_list(output_face)%cell%k
        face=output_mode_list(output_mode)%face_list(output_face)%point
 
! check the side of the surface on which we want output and change if required,
        if      (face.eq.face_xmin) then      
	  if (.NOT.output_mode_list(output_mode)%output_on_outward_normal) then ! output on other side of face
	    cx=cx-1
	    face=face_xmax
	  end if
         else if (face.eq.face_xmax) then	 
	  if (.NOT.output_mode_list(output_mode)%output_on_outward_normal) then ! output on other side of face
	    cx=cx+1
	    face=face_xmin
	  end if
        else if (face.eq.face_ymin) then	 
	  if (.NOT.output_mode_list(output_mode)%output_on_outward_normal) then ! output on other side of face
	    cy=cy-1
	    face=face_ymax
	  end if
        else if (face.eq.face_ymax) then	 
	  if (.NOT.output_mode_list(output_mode)%output_on_outward_normal) then ! output on other side of face
	    cy=cy+1
	    face=face_ymin
	  end if
        else if (face.eq.face_zmin) then	 
	  if (.NOT.output_mode_list(output_mode)%output_on_outward_normal) then ! output on other side of face
	    cz=cz-1
	    face=face_zmax
	  end if
        else if (face.eq.face_zmax) then	 
	  if (.NOT.output_mode_list(output_mode)%output_on_outward_normal) then ! output on other side of face
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

!! should the min face be copied to the output_surfaces(output_surface)%face_list? Use original face	

        output_mode_list(output_mode)%face_list(output_face)=output_face1
			 
      end do !next cell face in this surface	  
  
    end do ! next output_mode

  end if !n_output_modes.GT.0
  
  CALL write_line('FINISHED: set_mode_outputs_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_mode_outputs_in_mesh
!
! NAME
!     initialise_mode_outputs
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/01/2013 CJS
!
!
SUBROUTINE initialise_mode_outputs

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: output_mode
  integer	:: output_face
  integer	:: cx,cy,cz,face
   
! START
  
  CALL write_line('CALLED: initialise_mode_outputs',0,output_to_screen_flag)

  if (n_output_modes.GT.0) then
  
    do output_mode=1,n_output_modes
    
      write(*,*)'output_mode',output_mode
 
      do output_face=1,output_mode_list(output_mode)%n_mode_samples
      
        write(*,*)'output_face',output_face
      
! copy the faces from the geometry structure to the local
! output_surface structure
	
        cx=output_mode_list(output_mode)%face_list(output_face)%cell%i
        cy=output_mode_list(output_mode)%face_list(output_face)%cell%j
        cz=output_mode_list(output_mode)%face_list(output_face)%cell%k
        face=output_mode_list(output_mode)%face_list(output_face)%point

        write(*,*)'cell_face',cx,cy,cz,face
      
        if (rank.eq.cell_face_rank(cz,face)) then
   	     
          if      (face.eq.face_xmin) then      
            output_mode_list(output_mode)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_xmax) then	 
            output_mode_list(output_mode)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx+1,cy  ,cz  ,face_xmin)
          else if (face.eq.face_ymin) then	 
            output_mode_list(output_mode)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz  ,face_ymin)
          else if (face.eq.face_ymax) then	 
            output_mode_list(output_mode)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy+1,cz  ,face_ymin)
          else if (face.eq.face_zmin) then	 
            output_mode_list(output_mode)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz  ,face_zmin)
          else if (face.eq.face_zmax) then	 
            output_mode_list(output_mode)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz+1,face_zmin)
          end if
    					 
        end if ! output point belongs too this process
		 
      end do !next cell face in this surface	  
  
    end do ! next output_mode
  
    if (rank.eq.0) then
    
      OPEN(unit=mode_output_unit,file=trim(problem_name)//mode_output_extn)
  
      CALL write_time_domain_header_data(mode_output_unit,n_output_modes,n_timesteps)

    end if
    
  end if ! n_output_modes.GT.0
 
  CALL write_line('FINISHED: initialise_mode_outputs',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_mode_outputs
!
! SUBROUTINE Mode_output
!
! NAME
!     Mode_output
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/01/2013 CJS
!
!
SUBROUTINE Mode_output

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats

IMPLICIT NONE

! local variables

  integer	:: output_mode
  integer	:: output_face
  integer	:: cx,cy,cz,face,side

  integer 	:: number_of_faces
  integer 	:: field_component
  integer 	:: output_array_point
  real*8 	:: mode_field_value
  real*8 	:: value  
  real*8 	:: local_value  
  integer	:: n_data

! START
  
  CALL write_line('CALLED: Mode_output',0,timestepping_output_to_screen_flag)

  if (n_output_modes.GT.0) then
  
    do output_mode=1,n_output_modes
          
      number_of_faces=output_mode_list(output_mode)%n_mode_samples
      
      field_component=output_mode_list(output_mode)%field_component      
      
      value=0d0
      
      do output_face=1,output_mode_list(output_mode)%n_mode_samples
     
        cx  =output_mode_list(output_mode)%face_list(output_face)%cell%i
        cy  =output_mode_list(output_mode)%face_list(output_face)%cell%j
        cz  =output_mode_list(output_mode)%face_list(output_face)%cell%k
        face=output_mode_list(output_mode)%face_list(output_face)%point
      
        if (rank.eq.cell_face_rank(cz,face)) then
   
	  mode_field_value=output_mode_list(output_mode)%mode_field(output_face)
	  
          output_array_point=output_mode_list(output_mode)%face_output_field_number_list(output_face)
 	     
          if	(face.eq.face_xmin) then
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
	  
! Need to organise the correct side of the surface for the output field here...
	  value=value+face_output_field(output_array_point,side,field_component)*mode_field_value
	 	  
	end if ! output point in this processor's mesh
		  
      end do !next cell face in this surface	
 
! sum contributions from all processes to the rank 0 process    
#if defined(MPI)

      n_data=1
      CALL MPI_REDUCE(value, local_value,n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)

#elif defined(SEQ)

      local_value=value

#endif
  
      if (rank.eq.0) then
      
        if ( abs(local_value).lt.1D-30 )local_value=0d0

        write(mode_output_unit,time_domain_output_format)time,output_mode,local_value

      end if ! rank.eq.0
      
    end do ! next output surface

  end if ! n_output_surfaces.GT.0
  
  CALL write_line('FINISHED: Mode_output',0,timestepping_output_to_screen_flag)
  
  RETURN

END SUBROUTINE Mode_output
