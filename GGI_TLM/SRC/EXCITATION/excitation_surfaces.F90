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
! SUBROUTINE set_excitation_surfaces_in_mesh
!
! NAME
!     set_excitation_surfaces_in_mesh
!
! DESCRIPTION
!     loop over all the required excitations and flag all the excitation cells/ faces 
!     in the arrays local_cell_excitation(i,j,k) or local_surface_excitation(i,j,k,face)
!     as required
!     
! COMMENTS
!     Excitation_surfaces could be made more efficient in parallel - every process gets the whole excitation plane...
!  
!
! HISTORY
!
!     started 17/09/2012 CJS
!     Parallel 23/11/2012 CJS
!     Huygens surface 4/12/2012 CJS
!
!
SUBROUTINE set_excitation_surfaces_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information

IMPLICIT NONE

! local variables

  integer	:: excitation_point
  integer	:: excitation_surface
  integer	:: surface_number
  integer	:: number_of_faces
  integer	:: excitation_face
  integer	:: huygens_face
  integer	:: n_faces
  integer	:: cx,cy,cz,face
  type(cell_point)	:: excitation_face1
  type(cell_point)	:: excitation_face2
    
  logical inside

! START

  
  CALL write_line('CALLED: set_excitation_surfaces_in_mesh',0,output_to_screen_flag)
  
  if (n_excitation_surfaces.gt.0) then

! excitation SURFACES
    do excitation_surface=1,n_excitation_surfaces
    
      surface_number=excitation_surfaces(excitation_surface)%surface_number
      number_of_faces=problem_surfaces(surface_number)%number_of_faces
      excitation_surfaces(excitation_surface)%number_of_faces=number_of_faces
      
! allocate face list for the excitation surface      
      ALLOCATE( excitation_surfaces(excitation_surface)%face_list(1:number_of_faces) )
      
      do excitation_face=1,number_of_faces
      
! copy the faces from the geometry structure to the local
! excitation_surface structure
	
        cx=problem_surfaces(surface_number)%face_list(excitation_face)%cell%i
        cy=problem_surfaces(surface_number)%face_list(excitation_face)%cell%j
        cz=problem_surfaces(surface_number)%face_list(excitation_face)%cell%k
        face=problem_surfaces(surface_number)%face_list(excitation_face)%point
 
! check the side of the surface on which we want excitation and change if required,
        if      (face.eq.face_xmin) then      
	  if (.NOT.excitation_surfaces(excitation_surface)%excitation_on_outward_normal) then ! excitation on other side of face
	    cx=cx-1
	    face=face_xmax
	  end if
         else if (face.eq.face_xmax) then	 
	  if (.NOT.excitation_surfaces(excitation_surface)%excitation_on_outward_normal) then ! excitation on other side of face
	    cx=cx+1
	    face=face_xmin
	  end if
        else if (face.eq.face_ymin) then	 
	  if (.NOT.excitation_surfaces(excitation_surface)%excitation_on_outward_normal) then ! excitation on other side of face
	    cy=cy-1
	    face=face_ymax
	  end if
        else if (face.eq.face_ymax) then	 
	  if (.NOT.excitation_surfaces(excitation_surface)%excitation_on_outward_normal) then ! excitation on other side of face
	    cy=cy+1
	    face=face_ymin
	  end if
        else if (face.eq.face_zmin) then	 
	  if (.NOT.excitation_surfaces(excitation_surface)%excitation_on_outward_normal) then ! excitation on other side of face
	    cz=cz-1
	    face=face_zmax
	  end if
        else if (face.eq.face_zmax) then	 
	  if (.NOT.excitation_surfaces(excitation_surface)%excitation_on_outward_normal) then ! excitation on other side of face
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

        if (rank.eq.cell_face_rank(cz,face)) then
! excitation point belongs to this processor	  
          local_surface_excitation(cx  ,cy  ,cz  ,face)=1  
        end if ! excitation point belongs to this processor

!! should the min face be copied to the excitation_surfaces(excitation_surface)%face_list? Use original face	
!        excitation_surfaces(excitation_surface)%face_list(excitation_face)%cell%i=cx
!        excitation_surfaces(excitation_surface)%face_list(excitation_face)%cell%j=cy
!        excitation_surfaces(excitation_surface)%face_list(excitation_face)%cell%k=cz
!        excitation_surfaces(excitation_surface)%face_list(excitation_face)%point=face

         excitation_surfaces(excitation_surface)%face_list(excitation_face)=excitation_face1
			 
      end do !next cell face in this surface	  
  
    end do ! next excitation surface

  end if ! n_excitation_surfaces.GT.0
  
  CALL write_line('FINISHED: set_excitation_surfaces_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_excitation_surfaces_in_mesh
!
! NAME
!     initialise_excitation_surfaces
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
SUBROUTINE initialise_excitation_surfaces

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

  integer	:: excitation_surface
  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: excitation_face
  integer	:: huygens_face

  integer 	:: face_number
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4
    
! START

  
  CALL write_line('CALLED: initialise_excitation_surfaces',0,output_to_screen_flag)

! EXCITATION SURFACES

  if (n_excitation_surfaces.GT.0) then

    if (rank.eq.0) then
      write(info_file_unit,*)'Number of Excitation surfaces=',n_excitation_surfaces
    end if

    do excitation_surface=1,n_excitation_surfaces
    
      number_of_faces=excitation_surfaces(excitation_surface)%number_of_faces
      
      if (rank.eq.0) then
        write(info_file_unit,*)'Excitation surface number',excitation_surface,' Number of cell_faces=',number_of_faces
      end if
      
      ALLOCATE( excitation_surfaces(excitation_surface)%face_excitation_field_number_list(1:number_of_faces) )
      
      do excitation_face=1,number_of_faces
   
        cx  =excitation_surfaces(excitation_surface)%face_list(excitation_face)%cell%i
        cy  =excitation_surfaces(excitation_surface)%face_list(excitation_face)%cell%j
        cz  =excitation_surfaces(excitation_surface)%face_list(excitation_face)%cell%k
        face=excitation_surfaces(excitation_surface)%face_list(excitation_face)%point 
      
        if (rank.eq.cell_face_rank(cz,face)) then
   	     
          if      (face.eq.face_xmin) then      
            excitation_surfaces(excitation_surface)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_xmax) then	 
            excitation_surfaces(excitation_surface)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx+1,cy  ,cz  ,face_xmin)
          else if (face.eq.face_ymin) then	 
            excitation_surfaces(excitation_surface)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_ymin)
          else if (face.eq.face_ymax) then	 
            excitation_surfaces(excitation_surface)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy+1,cz  ,face_ymin)
          else if (face.eq.face_zmin) then	 
            excitation_surfaces(excitation_surface)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_zmin)
          else if (face.eq.face_zmax) then	 
            excitation_surfaces(excitation_surface)%face_excitation_field_number_list(excitation_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz+1,face_zmin)
          end if
    					 
        end if ! excitation point belongs too this process
		  
      end do !next cell face in this surface	  

    end do ! next excitation surface

  end if ! n_excitation_surfaces.GT.0 
  
  CALL write_line('FINISHED: initialise_excitation_surfaces',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_excitation_surfaces
!
! SUBROUTINE Surface_excitation
!
! NAME
!     Surface_excitation
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
SUBROUTINE Surface_excitation

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

  integer 	:: excitation_surface,excitation_face
  integer 	:: face,side
  integer 	:: cz
  integer 	:: number_of_faces
  integer 	:: function_number
  integer 	:: excitation_array_point
  integer 	:: field_component
  real*8	:: offset,offset_min
  real*8 	:: value  

! START
  
  CALL write_line('CALLED: Surface_excitation',0,timestepping_output_to_screen_flag)

  if (n_excitation_surfaces.gt.0) then

! excitation SURFACES
    do excitation_surface=1,n_excitation_surfaces
    
      function_number=excitation_surfaces(excitation_surface)%excitation_function_number
      field_component=excitation_surfaces(excitation_surface)%field_component
      value=excitation_functions(function_number)%value_face(timestep)
      
      number_of_faces=excitation_surfaces(excitation_surface)%number_of_faces
            
      do excitation_face=1,number_of_faces
      
        cz  =excitation_surfaces(excitation_surface)%face_list(excitation_face)%cell%k
        face  =excitation_surfaces(excitation_surface)%face_list(excitation_face)%point
      
        if (rank.eq.cell_face_rank(cz,face)) then
   
	  excitation_array_point=excitation_surfaces(excitation_surface)%face_excitation_field_number_list(excitation_face)
 	     
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
	
          face_excitation_field(excitation_array_point,side,field_component)=value
	  
	end if ! excitation point in this processor's mesh
		  
      end do !next cell face in this surface	  

    end do ! next excitation surface

  end if ! n_excitation_surfaces.GT.0
  
  CALL write_line('FINISHED: Surface_excitation',0,timestepping_output_to_screen_flag)
  
  RETURN

END SUBROUTINE Surface_excitation
