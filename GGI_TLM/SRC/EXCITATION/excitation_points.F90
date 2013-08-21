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
! SUBROUTINE set_excitation_points_in_mesh
!
! NAME
!     set_excitation_points_in_mesh
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
!
!
SUBROUTINE set_excitation_points_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information

IMPLICIT NONE

! local variables

  integer	:: excitation_point
  integer	:: cx,cy,cz,face
  type(cell_point)	:: excitation_face1
    
  logical inside

! START

  
  CALL write_line('CALLED: set_excitation_points_in_mesh',0,output_to_screen_flag)
  
  if (n_excitation_points.GT.0) then

    do excitation_point=1,n_excitation_points
  
      cx=excitation_points(excitation_point)%cell_point%cell%i
      cy=excitation_points(excitation_point)%cell_point%cell%j
      cz=excitation_points(excitation_point)%cell_point%cell%k
      face=excitation_points(excitation_point)%cell_point%point
  
      if (face.eq.centre) then
! cell centre excitation

        excitation_points(excitation_point)%rank =cell_rank(cz)

        if (rank.eq.cell_rank(cz)) then
! excitation point belongs to this processor
  
          local_cell_excitation(cx,cy,cz)=1
    
          write(*,*)'Setting cell centre excitation point',excitation_point
          write(*,*)'Coordinates:',cx,cy,cz
  
        end if ! output point belongs to this processor
  
      else
! must be an excitation point on a face
! set the excitation point number on the min face

        CALL get_min_face(excitation_points(excitation_point)%cell_point,excitation_face1)
	
        cx=excitation_face1%cell%i
        cy=excitation_face1%cell%j
        cz=excitation_face1%cell%k
	face=excitation_face1%point

        excitation_points(excitation_point)%rank =cell_rank(cz)

        if (rank.eq.cell_face_rank(cz,face)) then
! excitation point belongs to this processor
	  
          local_surface_excitation(cx  ,cy  ,cz  ,face)=1 
          write(*,*)'Setting cell face excitation point',excitation_point
          write(*,*)'Coordinates:',cx,cy,cz,' face:',face
  
        end if ! excitation point belongs to this processor
   
      end if ! centre or face excitation
  
    end do ! next excitation point

  end if ! n_excitation_points>0
  
  CALL write_line('FINISHED: set_excitation_points_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_excitation_points_in_mesh
!
! NAME
!     initialise_excitation_points
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
!     Huygens surface 4/12/2012 CJS
!
!
SUBROUTINE initialise_excitation_points

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

  integer	:: excitation_point
  integer	:: cx,cy,cz,face
    
! START
  
  CALL write_line('CALLED: initialise_excitation_points',0,output_to_screen_flag)

! EXCITATION_POINTS
  
  if (n_excitation_points.GT.0) then

    if (rank.eq.0) then
      write(info_file_unit,*)' ____________________________________________________'
      write(info_file_unit,*)''
      write(info_file_unit,*)'Number of Excitation points=',n_excitation_points
    end if
    
    do excitation_point=1,n_excitation_points
  
      cx=excitation_points(excitation_point)%cell_point%cell%i
      cy=excitation_points(excitation_point)%cell_point%cell%j
      cz=excitation_points(excitation_point)%cell_point%cell%k
      face=excitation_points(excitation_point)%cell_point%point
      
      if (rank.eq.cell_rank(cz)) then
  
        if (face.eq.centre) then
! cell centre excitation
  
          excitation_points(excitation_point)%cell_excitation_field_number=local_cell_excitation(cx,cy,cz)
  
        else
! must be an excitation point on a face
! set the excitation point number on the min face and give a -sign if it refers to the opposite side of the face
   	     
          if      (face.eq.face_xmin) then      
            excitation_points(excitation_point)%face_excitation_field_number=local_surface_excitation(cx  ,cy  ,cz  ,face_xmin)
          else if (face.eq.face_xmax) then	 
            excitation_points(excitation_point)%face_excitation_field_number=local_surface_excitation(cx+1,cy  ,cz  ,face_xmin)
          else if (face.eq.face_ymin) then	 
            excitation_points(excitation_point)%face_excitation_field_number=local_surface_excitation(cx  ,cy  ,cz  ,face_ymin)
          else if (face.eq.face_ymax) then	 
            excitation_points(excitation_point)%face_excitation_field_number=local_surface_excitation(cx  ,cy+1,cz  ,face_ymin)
          else if (face.eq.face_zmin) then	 
            excitation_points(excitation_point)%face_excitation_field_number=local_surface_excitation(cx  ,cy  ,cz  ,face_zmin)
          else if (face.eq.face_zmax) then	 
            excitation_points(excitation_point)%face_excitation_field_number=local_surface_excitation(cx  ,cy  ,cz+1,face_zmin)
          end if
  
        end if ! centre or face excitation
    					 
      end if ! excitation point belongs too this process
  
    end do ! next excitation point

  end if ! n_excitation_points>0
  
  CALL write_line('FINISHED: initialise_excitation_points',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_excitation_points
!
! SUBROUTINE Point_excitation
!
! NAME
!     Point_excitation
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
SUBROUTINE Point_excitation

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

  integer 	:: excitation_point
  integer 	:: face,side
  integer 	:: cz
  integer 	:: number_of_faces
  integer 	:: function_number
  integer 	:: excitation_array_point
  integer 	:: field_component
  real*8	:: offset,offset_min
  real*8 	:: value  

! START
  
  CALL write_line('CALLED: Point_excitation',0,timestepping_output_to_screen_flag)

! reset all fields in excitation array

  if (n_excitation_points.GT.0) then

    do excitation_point=1,n_excitation_points
    
      cz=excitation_points(excitation_point)%cell_point%cell%k
      face=excitation_points(excitation_point)%cell_point%point
      
      if (face.eq.centre) then
! cell centre excitation

        if (rank.eq.cell_rank(cz)) then
      
          function_number=excitation_points(excitation_point)%excitation_function_number
          field_component=excitation_points(excitation_point)%field_component

          excitation_array_point=excitation_points(excitation_point)%cell_excitation_field_number
          value=excitation_functions(function_number)%value(timestep)
          cell_excitation_field(excitation_array_point,field_component)=value
	  
	end if ! excitation point in this processor's mesh

      else
! must be an excitation point on a face

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

          excitation_array_point=excitation_points(excitation_point)%face_excitation_field_number
          value=excitation_functions(function_number)%value_face(timestep)
	
          face_excitation_field(excitation_array_point,side,field_component)=value
	
        end if ! excitation point in this processor's mesh
	
      end if ! centre or face excitation
  
    end do ! next excitation point

  end if ! n_excitation_points>0
  
  CALL write_line('FINISHED: Point_excitation',0,timestepping_output_to_screen_flag)
  
  RETURN

END SUBROUTINE Point_excitation
