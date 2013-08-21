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
!SUBROUTINE build_line_arc
!
! NAME
!     SUBROUTINE build_line_arc
!
! DESCRIPTION
!     build_line_arc:
!
!     create an arc
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 31/08/2012 CJS
!
!
SUBROUTINE build_line_arc(line_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: line_number

! local variables

real*8	:: radius,phi_min,phi_max
real*8	:: arc_length

integer	:: n_arc
integer :: n_phi
integer :: phi_loop
real*8	:: phi1,phi2,dphi

integer	:: number_of_line_segments
integer	:: line_segment_count

type(xyz)	:: point1,point2

! START

! CREATE AN ARC DEFINED BY RADIUS, PHI_MIN, PHI_MAX

      radius=problem_lines(line_number)%line_parameters(1)
      phi_min=problem_lines(line_number)%line_parameters(2)*pi/180d0
      phi_max=problem_lines(line_number)%line_parameters(3)*pi/180d0
      
      arc_length=radius*(phi_max-phi_min) 
          
      n_arc=2*NINT(arc_length/dl) 
      
      if (n_arc.lt.12) n_arc=12
      
      n_phi=n_arc
      dphi=(phi_max-phi_min)/n_phi
      
      number_of_line_segments=n_phi
      problem_lines(line_number)%number_of_line_segments=number_of_line_segments   
      allocate( problem_lines(line_number)%line_segment_list(1:number_of_line_segments) )

! create line segments    	  
      line_segment_count=0
      
      do phi_loop=1,n_phi
	
	phi1=phi_min+(phi_loop-1)*dphi
	phi2=phi_min+ phi_loop   *dphi
	  
! get points 1 and 2 on the arc
	CALL rphiz_to_xyz_point(radius,phi1,0d0,point1)
	CALL rphiz_to_xyz_point(radius,phi2,0d0,point2)

! apply the transformation to each of the points
        CALL apply_transformation(point1,problem_lines(line_number)%trans)
        CALL apply_transformation(point2,problem_lines(line_number)%trans)
	  	  
        line_segment_count=line_segment_count+1
        problem_lines(line_number)%line_segment_list(line_segment_count)%end(1)=point1
        problem_lines(line_number)%line_segment_list(line_segment_count)%end(2)=point2
		
      end do ! next phi

  RETURN

END SUBROUTINE build_line_arc
