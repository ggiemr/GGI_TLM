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
! SUBROUTINE set_mode_excitations_in_mesh
!
! NAME
!     set_mode_excitations_in_mesh
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
SUBROUTINE set_mode_excitations_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information

IMPLICIT NONE

! local variables

  integer	:: excitation_mode
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
  
  integer	:: excitation_face
  type(cell_point)	:: excitation_face1
  type(cell_point)	:: excitation_face2
  
! START
  
  CALL write_line('CALLED: set_mode_excitations_in_mesh',0,output_to_screen_flag)
  
  if (n_excitation_modes.GT.0) then
  
    do excitation_mode=1,n_excitation_modes
    
! get the coordinate range of the specified surface
      surface_number=excitation_mode_list(excitation_mode)%surface_number
      
      write(*,*)'Excitation mode=',excitation_mode
      write(*,*)'Surface number=',surface_number
     
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
      	   FILE=excitation_mode_list(excitation_mode)%mode_file_name,     &
      	   STATUS='old')
	   
! work out the number of field samples in the mode field file and the cell range in x, y and z
      n_field_samples=0
      
      mode_cxmin=nx
      mode_cxmax=0
      mode_cymin=ny
      mode_cymax=0
      mode_czmin=nz
      mode_czmax=0
      
      xcol=excitation_mode_list(excitation_mode)%xcol
      ycol=excitation_mode_list(excitation_mode)%ycol
      zcol=excitation_mode_list(excitation_mode)%zcol
      mode_col=excitation_mode_list(excitation_mode)%mode_col-3
      
      n_data_cols=excitation_mode_list(excitation_mode)%mode_col-3
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
	
	write(*,*)'Error in set_mode_excitations_in_mesh:'
	write(*,*)'Mode surface is not a plane'
	write(*,*)' mode_nx=',mode_nx
	write(*,*)' mode_ny=',mode_ny
	write(*,*)' mode_nz=',mode_nz
	STOP
	   
      end if
      
! check that excitation surface is a plane

      if ( (surface_nx.ne.1).AND.	&
           (surface_ny.ne.1).AND.	&
           (surface_nz.ne.1) ) then
	
	write(*,*)'Error in set_mode_excitations_in_mesh:'
	write(*,*)'mode excitation surface is not a plane'
	write(*,*)' surface_nx=',surface_nx
	write(*,*)' surface_ny=',surface_ny
	write(*,*)' surface_nz=',surface_nz
	STOP
	   
      end if
      
! check that the cell ranges for surface and mode are the same

      if ( (mode_nx.ne.surface_nx).OR.	&
           (mode_ny.ne.surface_ny).OR.	&
           (mode_nz.ne.surface_nz) ) then
	
	write(*,*)'Error in set_mode_excitations_in_mesh:'
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
      
! calculate the number of mode excitation cells
      n_surface_cells=surface_nx*surface_ny*surface_nz
      n_mode_cells=mode_nx*mode_ny*mode_nz

! check consistency      
      if (n_surface_cells.ne.n_mode_cells) then
	
	write(*,*)'Error in set_mode_excitations_in_mesh:'
	write(*,*)'Discrepancy between number of mode surface cells and excitation surface cells'
	write(*,*)'n_surface_cells=',n_surface_cells,' n_mode_cells=',n_mode_cells
	STOP
	   
      end if
      
! Removed check - fails in parallel as problem_surfaces(surface_number)%number_of_faces
! is the number of faces belonging to this process only.
!
! check the number of surface cells
!      if (n_mode_cells.ne.problem_surfaces(surface_number)%number_of_faces) then
!	
!	write(*,*)'Error in set_mode_excitations_in_mesh:'
!	write(*,*)'Discrepancy in number of excitation surface cells'
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
	
	    write(*,*)'Error in set_mode_excitations_in_mesh:'
	    write(*,*)'Not all the surface cell faces are the same'
	    write(*,*)'face(1)=',face
	    write(*,*)'i=',i
	    write(*,*)'face(i)=',problem_surfaces(surface_number)%face_list(i)%point
	    STOP
	  
	  end if
        end do ! next surface cell face

        if ( (face.eq.face_xmin).OR.(face.eq.face_xmax) ) then
          if (mode_nx.ne.1) then
	    write(*,*)'Error in set_mode_excitations_in_mesh:'
	    write(*,*)'face normal=x and mode_nx.ne.1'
	    STOP
	  end if
        else if ( (face.eq.face_ymin).OR.(face.eq.face_ymax) ) then
          if (mode_ny.ne.1) then
	    write(*,*)'Error in set_mode_excitations_in_mesh:'
	    write(*,*)'face normal=y and mode_ny.ne.1'
	    STOP
	  end if     
        else if ( (face.eq.face_zmin).OR.(face.eq.face_zmax) ) then
          if (mode_nz.ne.1) then
	    write(*,*)'Error in set_mode_excitations_in_mesh:'
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

! allocate memory for the mode_excitation      

        excitation_mode_list(excitation_mode)%n_mode_samples=n_mode_samples
  
        ALLOCATE( excitation_mode_list(excitation_mode)%face_list(1:n_mode_samples) )
        ALLOCATE( excitation_mode_list(excitation_mode)%mode_field(1:n_mode_samples) )
        ALLOCATE( excitation_mode_list(excitation_mode)%face_excitation_field_number_list(1:n_mode_samples) )
      
! re-read the mode file data into the mode_excitation structure
     
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
	    excitation_mode_list(excitation_mode)%face_list(n_mode_samples)%cell%i=cx
	    excitation_mode_list(excitation_mode)%face_list(n_mode_samples)%cell%j=cy
	    excitation_mode_list(excitation_mode)%face_list(n_mode_samples)%cell%k=cz
	    excitation_mode_list(excitation_mode)%face_list(n_mode_samples)%point=face
	    excitation_mode_list(excitation_mode)%mode_field(n_mode_samples)=Field
	  end if
	
        end do ! next sample
      
      else
! no face samples in this process

        excitation_mode_list(excitation_mode)%n_mode_samples= 0  
	   
      end if ! n_faces.ne.0     
	       
      DEALLOCATE( data_cols )
    
      close(unit=mode_file_unit)

! set the mode excitation faces in the mesh     
 
      do excitation_face=1,excitation_mode_list(excitation_mode)%n_mode_samples
      
! copy the faces from the geometry structure to the local
! excitation_surface structure
	
        cx=excitation_mode_list(excitation_mode)%face_list(excitation_face)%cell%i
        cy=excitation_mode_list(excitation_mode)%face_list(excitation_face)%cell%j
        cz=excitation_mode_list(excitation_mode)%face_list(excitation_face)%cell%k
        face=excitation_mode_list(excitation_mode)%face_list(excitation_face)%point
 
! check the side of the surface on which we want excitation and change if required,
        if      (face.eq.face_xmin) then      
	  if (.NOT.excitation_mode_list(excitation_mode)%excitation_on_outward_normal) then ! excitation on other side of face
	    cx=cx-1
	    face=face_xmax
	  end if
         else if (face.eq.face_xmax) then	 
	  if (.NOT.excitation_mode_list(excitation_mode)%excitation_on_outward_normal) then ! excitation on other side of face
	    cx=cx+1
	    face=face_xmin
	  end if
        else if (face.eq.face_ymin) then	 
	  if (.NOT.excitation_mode_list(excitation_mode)%excitation_on_outward_normal) then ! excitation on other side of face
	    cy=cy-1
	    face=face_ymax
	  end if
        else if (face.eq.face_ymax) then	 
	  if (.NOT.excitation_mode_list(excitation_mode)%excitation_on_outward_normal) then ! excitation on other side of face
	    cy=cy+1
	    face=face_ymin
	  end if
        else if (face.eq.face_zmin) then	 
	  if (.NOT.excitation_mode_list(excitation_mode)%excitation_on_outward_normal) then ! excitation on other side of face
	    cz=cz-1
	    face=face_zmax
	  end if
        else if (face.eq.face_zmax) then	 
	  if (.NOT.excitation_mode_list(excitation_mode)%excitation_on_outward_normal) then ! excitation on other side of face
	    cz=cz+1
	    face=face_zmin
	  end if
        end if
  	     
! Set the excitation face in the local_surface_excitation array
! We must set the excitation point number on the min face 

        excitation_face1%cell%i=cx
        excitation_face1%cell%j=cy
        excitation_face1%cell%k=cz
        excitation_face1%point=face
	
        CALL get_min_face(excitation_face1,excitation_face2)
	
        cx=excitation_face2%cell%i
        cy=excitation_face2%cell%j
        cz=excitation_face2%cell%k
	face=excitation_face2%point
	
        if (rank.eq.cell_rank(cz)) then
! excitation point belongs to this processor	  
          local_surface_excitation(cx  ,cy  ,cz  ,face)=1  
        end if ! excitation point belongs to this processor

!! should the min face be copied to the excitation_surfaces(excitation_surface)%face_list? Use original face	

        excitation_mode_list(excitation_mode)%face_list(excitation_face)=excitation_face1
			 
      end do !next cell face in this surface	  
  
    end do ! next excitation_mode

  end if !n_excitation_modes.GT.0
  
  CALL write_line('FINISHED: set_mode_excitations_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_mode_excitations_in_mesh
!
! NAME
!     initialise_mode_excitations
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
SUBROUTINE initialise_mode_excitations

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: excitation_mode
  integer	:: excitation_face
  integer	:: cx,cy,cz,face
   
! START
  
  CALL write_line('CALLED: initialise_mode_excitations',0,output_to_screen_flag)

  if (n_excitation_modes.GT.0) then
  
    do excitation_mode=1,n_excitation_modes
 
      do excitation_face=1,excitation_mode_list(excitation_mode)%n_mode_samples
      
! copy the faces from the geometry structure to the local
! excitation_surface structure
	
        cx=excitation_mode_list(excitation_mode)%face_list(excitation_face)%cell%i
        cy=excitation_mode_list(excitation_mode)%face_list(excitation_face)%cell%j
        cz=excitation_mode_list(excitation_mode)%face_list(excitation_face)%cell%k
        face=excitation_mode_list(excitation_mode)%face_list(excitation_face)%point
      
        if (rank.eq.cell_face_rank(cz,face)) then
   	     
          if      (face.eq.face_xmin) then      
            excitation_mode_list(excitation_mode)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_xmax) then	 
            excitation_mode_list(excitation_mode)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx+1,cy  ,cz  ,face_xmin)
          else if (face.eq.face_ymin) then	 
            excitation_mode_list(excitation_mode)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_ymin)
          else if (face.eq.face_ymax) then	 
            excitation_mode_list(excitation_mode)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy+1,cz  ,face_ymin)
          else if (face.eq.face_zmin) then	 
            excitation_mode_list(excitation_mode)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_zmin)
          else if (face.eq.face_zmax) then	 
            excitation_mode_list(excitation_mode)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz+1,face_zmin)
          end if
    					 
        end if ! excitation point belongs too this process
			 
      end do !next cell face in this surface	  
  
    end do ! next excitation_mode

  end if ! n_excitation_modes.GT.0
 
  CALL write_line('FINISHED: initialise_mode_excitations',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_mode_excitations
!
! SUBROUTINE Mode_excitation
!
! NAME
!     Mode_excitation
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
SUBROUTINE Mode_excitation

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

  integer	:: excitation_mode,excitation_face
  integer	:: cz,face,side

  integer 	:: number_of_faces
  integer 	:: function_number
  integer 	:: excitation_array_point
  integer 	:: field_component
  real*8 	:: mode_field_value
  real*8 	:: value  

! START
  
  CALL write_line('CALLED: Mode_excitation',0,timestepping_output_to_screen_flag)

  if (n_excitation_modes.GT.0) then
  
    do excitation_mode=1,n_excitation_modes
    
      function_number=excitation_mode_list(excitation_mode)%excitation_function_number
      field_component=excitation_mode_list(excitation_mode)%field_component
      value=excitation_functions(function_number)%value_face(timestep)
      
      number_of_faces=excitation_mode_list(excitation_mode)%n_mode_samples
            
      do excitation_face=1,excitation_mode_list(excitation_mode)%n_mode_samples
     
        cz  =excitation_mode_list(excitation_mode)%face_list(excitation_face)%cell%k
        face  =excitation_mode_list(excitation_mode)%face_list(excitation_face)%point
      
        if (rank.eq.cell_face_rank(cz,face)) then
 	     
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
   
	  excitation_array_point=excitation_mode_list(excitation_mode)%face_excitation_field_number_list(excitation_face)
	  mode_field_value=excitation_mode_list(excitation_mode)%mode_field(excitation_face)
	  
          face_excitation_field(excitation_array_point,side,field_component)=mode_field_value*value
	  	  
	end if ! excitation point in this processor's mesh
		  
      end do !next cell face in this surface	  

    end do ! next excitation surface

  end if ! n_excitation_surfaces.GT.0
  
  CALL write_line('FINISHED: Mode_excitation',0,timestepping_output_to_screen_flag)
  
  RETURN

END SUBROUTINE Mode_excitation
