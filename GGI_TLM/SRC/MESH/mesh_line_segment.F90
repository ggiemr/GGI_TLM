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
!SUBROUTINE mesh_line_segment
!
! NAME
!     SUBROUTINE mesh_line_segment
!
! DESCRIPTION
!     mesh_line_segment:
!
!     mesh a single line_segment. The resulting segments run from cell centre to cell centre. 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 31/08/2012 CJS
!
!
  SUBROUTINE mesh_line_segment(line_number,segment_number,segment_count,set_mesh_flag)
 
USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE cell_parameters
USE constants

IMPLICIT NONE

integer				:: line_number
integer				:: segment_number
integer				:: segment_count
logical				:: set_mesh_flag

! local variables

  type(ijk)		:: cell1,cell2
  
  integer 		:: ix,iy,iz
  
  integer 		:: n_segments
  integer 		:: n_cells

! START
  
! get line segment end point cells
  
  call point_to_cell(problem_lines(line_number)%line_segment_list(segment_number)%end(1),cell1)
  call point_to_cell(problem_lines(line_number)%line_segment_list(segment_number)%end(2),cell2)
  
  if ( (cell1%i.EQ.cell2%i).AND. (cell1%j.EQ.cell2%j).AND. (cell1%k.EQ.cell2%k) ) then
! both ends of the line sement mesh to the same point so return

    RETURN  
    
  end if
  
  if ( (abs(cell1%i-cell2%i).GT.1).OR.	&
       (abs(cell1%j-cell2%j).GT.1).OR.	&
       (abs(cell1%k-cell2%k).GT.1) ) then
! line segment is longer than cell dimension, this causes an error for the mesh generation

    GOTO 9000
    
  end if
  
  n_segments=0
  
  if ( (cell2%i-cell1%i).EQ.1 ) then
! step in + x direction from cell 1 to cell 2
    
    if(set_mesh_flag) then      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%i=cell1%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%point=centre
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%i=cell1%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%point=face_xmax
      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%point=face_xmin
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%point=centre
    end if ! set_mesh_flag

    n_segments=n_segments+2

  end if
  
  if ( (cell2%i-cell1%i).EQ.-1 ) then
! step in - x direction from cell 1 to cell 2
    
    if(set_mesh_flag) then      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%i=cell1%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%point=centre
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%i=cell1%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%point=face_xmin
      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%point=face_xmax
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%point=centre
    end if ! set_mesh_flag

    n_segments=n_segments+2

  end if
  
  if ( (cell2%j-cell1%j).EQ.1 ) then
! step in + y direction from cell 1 to cell 2
    
    if(set_mesh_flag) then

      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%point=centre
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%point=face_ymax
      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%point=face_ymin
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%point=centre

    end if ! set_mesh_flag

    n_segments=n_segments+2

  end if
  
  if ( (cell2%j-cell1%j).EQ.-1 ) then
! step in - y direction from cell 1 to cell 2
    
    if(set_mesh_flag) then

      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%point=centre
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%j=cell1%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%point=face_ymin
      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%point=face_ymax
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%point=centre

    end if ! set_mesh_flag

    n_segments=n_segments+2

  end if
  
  if ( (cell2%k-cell1%k).EQ.1 ) then
! step in + z direction from cell 1 to cell 2
    
    if(set_mesh_flag) then

      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%point=centre
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%point=face_zmax
      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%k=cell2%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%point=face_zmin
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%k=cell2%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%point=centre

    end if ! set_mesh_flag

    n_segments=n_segments+2

  end if
  
  if ( (cell2%k-cell1%k).EQ.-1 ) then
! step in - z direction from cell 1 to cell 2
    
    if(set_mesh_flag) then

      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(1)%point=centre
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%cell%k=cell1%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+1)%segment_point(2)%point=face_zmin
      
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%cell%k=cell2%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(1)%point=face_zmax
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%i=cell2%i
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%j=cell2%j
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%cell%k=cell2%k
      problem_lines(line_number)%cell_segment_list(segment_count+n_segments+2)%segment_point(2)%point=centre
      
    end if ! set_mesh_flag

    n_segments=n_segments+2

  end if
  
  segment_count=segment_count+n_segments
  
!  write(*,*)'FINISHED: mesh_line_segment',n_cells,segment_count

  RETURN

9000 CALL write_line('Error in mesh_line_segment',0,.TRUE.)
     CALL write_line('Line segment ',0,.TRUE.)
  
  END SUBROUTINE mesh_line_segment
