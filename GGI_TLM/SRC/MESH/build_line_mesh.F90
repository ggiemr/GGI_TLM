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
!SUBROUTINE build_line_mesh
!
! NAME
!     SUBROUTINE build_line_mesh
!
! DESCRIPTION
!     build_line_mesh:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE build_line_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: line_number

integer	:: number_of_line_segments

integer		:: segment_number
integer		:: segment_count
logical		:: set_mesh_flag

! START

  CALL write_line('CALLED: build_line_mesh',0,output_to_screen_flag)

  do line_number=1,n_lines

    number_of_line_segments=problem_lines(line_number)%number_of_line_segments

! Initially only count the number of cell segments in the line mesh

    set_mesh_flag=.FALSE.
    segment_count=0
    do segment_number=1,number_of_line_segments  
      CALL mesh_line_segment(line_number,segment_number,segment_count,set_mesh_flag)   
    end do ! next line segment

! Allocate memory for the mesh

    problem_lines(line_number)%number_of_cell_segments=segment_count
    
    if (segment_count.gt.0) then
    
      ALLOCATE( problem_lines(line_number)%cell_segment_list(1:segment_count) )
    
! Now regenerate the line mesh and set the cell_segment_list

      set_mesh_flag=.TRUE.
      segment_count=0
      do segment_number=1,number_of_line_segments  
        CALL mesh_line_segment(line_number,segment_number,segment_count,set_mesh_flag)   
      end do ! next line segment
      
    end if ! number of segments .gt.0
    
  end do ! next line number

  CALL write_line('FINISHED: build_line_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE build_line_mesh
