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
!SUBROUTINE trim_mesh
!
! NAME
!     SUBROUTINE trim_mesh
!
! DESCRIPTION
!     trim_mesh:
!
!     once we have read volume, surface, line and point meshes from file, remove all
!     elements which are outside the defined problem boundary
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 12/09/2012 CJS
!
!
SUBROUTINE trim_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE mesh
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number
integer	:: cell,number_of_cells,new_number_of_cells

integer	:: surface_number
integer	:: face,number_of_faces,new_number_of_faces

integer	:: line_number
integer	:: segment,number_of_cell_segments,new_number_of_cell_segments

integer	:: point_number,new_number_of_points

type(cell_point)	:: local_cell_point
logical			:: inside
logical			:: inside1,inside2

type(cell_point),allocatable	:: new_cell_point_list(:)
type(cell_segment),allocatable	:: new_cell_segment_list(:)

! START

  CALL write_line('CALLED: trim_mesh',0,output_to_screen_flag)
  
! TRIM VOLUMES
  do volume_number=1,n_volumes
    
    number_of_cells=problem_volumes(volume_number)%number_of_cells
    
    ALLOCATE( new_cell_point_list(1:number_of_cells) )
    
    new_number_of_cells=0
    do cell=1,number_of_cells
    
      CALL cell_point_inside_mesh(problem_volumes(volume_number)%cell_list(cell), inside)
				  
      if (inside) then
      
        new_number_of_cells=new_number_of_cells+1
	
	new_cell_point_list(new_number_of_cells)%cell%i=problem_volumes(volume_number)%cell_list(cell)%cell%i
	new_cell_point_list(new_number_of_cells)%cell%j=problem_volumes(volume_number)%cell_list(cell)%cell%j
	new_cell_point_list(new_number_of_cells)%cell%k=problem_volumes(volume_number)%cell_list(cell)%cell%k
	new_cell_point_list(new_number_of_cells)%point =problem_volumes(volume_number)%cell_list(cell)%point
	
      end if ! cell_point inside mesh
      
    end do ! next cell
    
    CALL write_line_integer('Volume number',volume_number,0,output_to_screen_flag)
    CALL write_line_integer('Old number of cells',number_of_cells,0,output_to_screen_flag)
    CALL write_line_integer('New number of cells',new_number_of_cells,0,output_to_screen_flag)

! reallocate cell list and fill with the OK cells

    DEALLOCATE( problem_volumes(volume_number)%cell_list )
    
    number_of_cells=new_number_of_cells
    problem_volumes(volume_number)%number_of_cells=number_of_cells
    
    ALLOCATE( problem_volumes(volume_number)%cell_list(1:number_of_cells) )
    
    do cell=1,number_of_cells
      problem_volumes(volume_number)%cell_list(cell)%cell%i=new_cell_point_list(new_number_of_cells)%cell%i
      problem_volumes(volume_number)%cell_list(cell)%cell%j=new_cell_point_list(new_number_of_cells)%cell%j
      problem_volumes(volume_number)%cell_list(cell)%cell%k=new_cell_point_list(new_number_of_cells)%cell%k
      problem_volumes(volume_number)%cell_list(cell)%point =new_cell_point_list(new_number_of_cells)%point
    end do ! next cell
    
    DEALLOCATE( new_cell_point_list )
    
  end do ! next volume number
 
! TRIM SURFACES
  do surface_number=1,n_surfaces
    
    number_of_faces=problem_surfaces(surface_number)%number_of_faces
       
    ALLOCATE( new_cell_point_list(1:number_of_faces) )
    
    new_number_of_faces=0
    
    do face=1,number_of_faces
    
      CALL cell_point_inside_mesh(problem_surfaces(surface_number)%face_list(face),inside)
      if (inside) then
      
        new_number_of_faces=new_number_of_faces+1
	
	new_cell_point_list(new_number_of_faces)%cell%i=problem_surfaces(surface_number)%face_list(face)%cell%i
	new_cell_point_list(new_number_of_faces)%cell%j=problem_surfaces(surface_number)%face_list(face)%cell%j
	new_cell_point_list(new_number_of_faces)%cell%k=problem_surfaces(surface_number)%face_list(face)%cell%k
	new_cell_point_list(new_number_of_faces)%point =problem_surfaces(surface_number)%face_list(face)%point
	
      end if ! cell_point inside mesh

    end do ! next face
    
    CALL write_line_integer('Surface number',surface_number,0,output_to_screen_flag)
    CALL write_line_integer('Old number of faces',number_of_faces,0,output_to_screen_flag)
    CALL write_line_integer('New number of faces',new_number_of_faces,0,output_to_screen_flag)

! reallocate face list and fill with the OK cells

    DEALLOCATE( problem_surfaces(surface_number)%face_list )
    
    number_of_faces=new_number_of_faces
    problem_surfaces(surface_number)%number_of_faces=number_of_faces
    
    ALLOCATE( problem_surfaces(surface_number)%face_list(1:number_of_faces) )
    
    do face=1,number_of_faces
    
      problem_surfaces(surface_number)%face_list(face)%cell%i=new_cell_point_list(face)%cell%i
      problem_surfaces(surface_number)%face_list(face)%cell%j=new_cell_point_list(face)%cell%j  
      problem_surfaces(surface_number)%face_list(face)%cell%k=new_cell_point_list(face)%cell%k  
      problem_surfaces(surface_number)%face_list(face)%point =new_cell_point_list(face)%point	
      
    end do ! next cell
    
    DEALLOCATE( new_cell_point_list )
    
  end do ! next surface number
 
! TRIM LINES
  do line_number=1,n_lines

    number_of_cell_segments=problem_lines(line_number)%number_of_cell_segments
       
    ALLOCATE( new_cell_segment_list(1:number_of_cell_segments) )
    
    new_number_of_cell_segments=0
        
    do segment=1,number_of_cell_segments
    
      CALL cell_point_inside_mesh(problem_lines(line_number)%cell_segment_list(segment)%segment_point(1),inside1 )
      CALL cell_point_inside_mesh(problem_lines(line_number)%cell_segment_list(segment)%segment_point(2),inside2 )
				
      if (inside1.AND.inside2) then
      
        new_number_of_cell_segments=new_number_of_cell_segments+1
	
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(1)%cell%i=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%i
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(1)%cell%j=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%j
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(1)%cell%k=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%k
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(1)%point=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%point
	
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(2)%cell%i=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%i
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(2)%cell%j=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%j
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(2)%cell%k=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%k
	new_cell_segment_list(new_number_of_cell_segments)%segment_point(2)%point=	&
	  problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%point	
	
      end if ! cell_point inside mesh
      
    end do ! next segment
    
    CALL write_line_integer('Line number',line_number,0,output_to_screen_flag)
    CALL write_line_integer('Old number of cell segments',number_of_cell_segments,0,output_to_screen_flag)
    CALL write_line_integer('New number of cell segments',new_number_of_cell_segments,0,output_to_screen_flag)
      
! reallocate cell segment list and fill with the OK cell segments

    DEALLOCATE( problem_lines(line_number)%cell_segment_list )
    
    number_of_cell_segments=new_number_of_cell_segments
    problem_lines(line_number)%number_of_cell_segments=number_of_cell_segments
    
    ALLOCATE( problem_lines(line_number)%cell_segment_list(1:number_of_cell_segments) )
    
    do segment=1,number_of_cell_segments
    
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%i=	&
	                       new_cell_segment_list(segment)%segment_point(1)%cell%i
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%j=	&
	                       new_cell_segment_list(segment)%segment_point(1)%cell%j
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%k=	&
	                       new_cell_segment_list(segment)%segment_point(1)%cell%k
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%point =	&
	                       new_cell_segment_list(segment)%segment_point(1)%point
    
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%i=	&
	                       new_cell_segment_list(segment)%segment_point(2)%cell%i
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%j=	&
	                       new_cell_segment_list(segment)%segment_point(2)%cell%j
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%k=	&
	                       new_cell_segment_list(segment)%segment_point(2)%cell%k
 	problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%point =	&
	                       new_cell_segment_list(segment)%segment_point(2)%point
				
    end do ! next cell segment
       
    DEALLOCATE( new_cell_segment_list )
    
  end do ! next line number

! TRIM POINTS  

  new_number_of_points=0
  
  ALLOCATE( new_cell_point_list(1:n_points) )
 
  do point_number=1,n_points

    local_cell_point%cell%i=problem_points(point_number)%cell%i
    local_cell_point%cell%j=problem_points(point_number)%cell%j
    local_cell_point%cell%k=problem_points(point_number)%cell%k
    local_cell_point%point =centre
    
    CALL cell_point_inside_mesh(local_cell_point,inside )
    
    if (inside) then
    
      new_number_of_points=new_number_of_points+1

      new_cell_point_list(new_number_of_points)%cell%i=local_cell_point%cell%i
      new_cell_point_list(new_number_of_points)%cell%j=local_cell_point%cell%j
      new_cell_point_list(new_number_of_points)%cell%k=local_cell_point%cell%k

    end if ! cell_point inside mesh
    
  end do ! next point number
    
  CALL write_line('Point list',0,output_to_screen_flag)
  CALL write_line_integer('Old number of points',n_points,0,output_to_screen_flag)
  CALL write_line_integer('New number of points',new_number_of_points,0,output_to_screen_flag)
 
! reallocate point list  
  n_points=new_number_of_points
  
  DEALLOCATE( problem_points )
  
  ALLOCATE( problem_points(1:n_points) )
  
  do point_number=1,n_points

    problem_points(point_number)%cell%i=new_cell_point_list(point_number)%cell%i
    problem_points(point_number)%cell%j=new_cell_point_list(point_number)%cell%j
    problem_points(point_number)%cell%k=new_cell_point_list(point_number)%cell%k
        
  end do ! next point number
  
  DEALLOCATE( new_cell_point_list )
  
  CALL write_line('FINISHED: trim_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE trim_mesh
