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
!SUBROUTINE build_line_rectangle
!
! NAME
!     SUBROUTINE build_line_rectangle
!
! DESCRIPTION
!     build_line_rectangle:
!
!     create a straight line defined by its end points
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
SUBROUTINE build_line_rectangle(line_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: line_number

! local variables

real*8  :: lx2,ly2
real*8	:: line_length(5)
real*8	:: xmin(5),ymin(5),zmin(5),xmax(5),ymax(5),zmax(5),dx(5),dy(5),dz(5)

integer	:: line_segment_loop

integer	:: number_of_line_segments(5),total_number_of_line_segments
integer	:: line_segment_count

type(xyz)	:: point1,point2

integer :: line

! START

! CREATE A LINE SEGMENT DEFINED BY LINE LENGTH AND TRANSFORMATION

    lx2=problem_lines(line_number)%line_parameters(1)/2d0
    ly2=problem_lines(line_number)%line_parameters(2)/2d0
    
    zmin(:)=0d0
    zmax(:)=0d0
    
    total_number_of_line_segments=0

    do line=1,5
    
      if (line.EQ.1) then
      
        xmin(line)= lx2
        xmax(line)= lx2
        ymin(line)= 0d0
        ymax(line)= ly2

      else if (line.EQ.2) then
      
        xmin(line)= lx2
        xmax(line)=-lx2
        ymin(line)= ly2
        ymax(line)= ly2
    
      else if (line.EQ.3) then
      
        xmin(line)=-lx2
        xmax(line)=-lx2
        ymin(line)= ly2
        ymax(line)=-ly2

      else if (line.EQ.4) then
      
        xmin(line)=-lx2
        xmax(line)= lx2
        ymin(line)=-ly2
        ymax(line)=-ly2

      else if (line.EQ.5) then
      
        xmin(line)= lx2
        xmax(line)= lx2
        ymin(line)=-ly2
        ymax(line)= 0d0
        
      end if
      
      line_length(line)=sqrt( (xmax(line)-xmin(line))**2+(ymax(line)-ymin(line))**2+(zmax(line)-zmin(line))**2 ) 
      
      number_of_line_segments(line)=2*INT(line_length(line)/dl)
      if (number_of_line_segments(line).lt.4) number_of_line_segments(line)=4 
      
      total_number_of_line_segments=total_number_of_line_segments+number_of_line_segments(line)
   
      dx(line)=(xmax(line)-xmin(line))/number_of_line_segments(line)
      dy(line)=(ymax(line)-ymin(line))/number_of_line_segments(line)
      dz(line)=(zmax(line)-zmin(line))/number_of_line_segments(line)
      
    end do  ! next line
        
    problem_lines(line_number)%number_of_line_segments=total_number_of_line_segments   
    allocate( problem_lines(line_number)%line_segment_list(1:total_number_of_line_segments) )
	  
! create line segments    	  
    line_segment_count=0
    
    do line=1,5
      
      do line_segment_loop=1,number_of_line_segments(line)

        point1%x=xmin(line)+(line_segment_loop-1)*dx(line)
        point1%y=ymin(line)+(line_segment_loop-1)*dy(line)
        point1%z=zmin(line)+(line_segment_loop-1)*dz(line)
        point2%x=xmin(line)+line_segment_loop*dx(line)
        point2%y=ymin(line)+line_segment_loop*dy(line)
        point2%z=zmin(line)+line_segment_loop*dz(line)
	
	if (line_segment_loop.EQ.number_of_line_segments(line)) then
          point2%x=xmax(line)
          point2%y=ymax(line)
          point2%z=zmax(line)
	end if
    
! apply the transformation to each of the points
        CALL apply_transformation(point1,problem_lines(line_number)%trans)
        CALL apply_transformation(point2,problem_lines(line_number)%trans)
	  
        line_segment_count=line_segment_count+1
        problem_lines(line_number)%line_segment_list(line_segment_count)%end(1)=point1
        problem_lines(line_number)%line_segment_list(line_segment_count)%end(2)=point2
	
      end do ! next line segment
	
    end do ! next line

  RETURN

END SUBROUTINE build_line_rectangle
