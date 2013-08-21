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
!SUBROUTINE build_line_straight_line
!
! NAME
!     SUBROUTINE build_line_straight_line
!
! DESCRIPTION
!     build_line_straight_line:
!
!     create a straight line
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
SUBROUTINE build_line_straight_line(line_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: line_number

! local variables

real*8	:: line_length
real*8	:: xmin,ymin,zmin,xmax,ymax,zmax,dx,dy,dz

integer	:: line_segment_loop

integer	:: number_of_line_segments
integer	:: line_segment_count

type(xyz)	:: point1,point2

! START

! CREATE A LINE SEGMENT DEFINED BY LINE LENGTH AND TRANSFORMATION

      xmin=-problem_lines(line_number)%line_parameters(1)/2d0
      ymin=0d0
      zmin=0d0
      xmax=problem_lines(line_number)%line_parameters(1)/2d0
      ymax=0d0
      zmax=0d0
      
      line_length=sqrt( (xmax-xmin)**2+(ymax-ymin)**2+(zmax-zmin)**2 ) 
      
      number_of_line_segments=2*INT(line_length/dl)
      if (number_of_line_segments.lt.4) number_of_line_segments=4 
   
      dx=(xmax-xmin)/number_of_line_segments
      dy=(ymax-ymin)/number_of_line_segments
      dz=(zmax-zmin)/number_of_line_segments
   
      problem_lines(line_number)%number_of_line_segments=number_of_line_segments   
      allocate( problem_lines(line_number)%line_segment_list(1:number_of_line_segments) )
	  
! create line segments    	  
      line_segment_count=0
      
      do line_segment_loop=1,number_of_line_segments

        point1%x=xmin+(line_segment_loop-1)*dx
        point1%y=ymin+(line_segment_loop-1)*dy
        point1%z=zmin+(line_segment_loop-1)*dz
        point2%x=xmin+line_segment_loop*dx
        point2%y=ymin+line_segment_loop*dy
        point2%z=zmin+line_segment_loop*dz
    
! apply the transformation to each of the points
        CALL apply_transformation(point1,problem_lines(line_number)%trans)
        CALL apply_transformation(point2,problem_lines(line_number)%trans)
	  
        line_segment_count=line_segment_count+1
        problem_lines(line_number)%line_segment_list(line_segment_count)%end(1)=point1
        problem_lines(line_number)%line_segment_list(line_segment_count)%end(2)=point2
	
      end do ! next line segment
      
  RETURN

END SUBROUTINE build_line_straight_line
