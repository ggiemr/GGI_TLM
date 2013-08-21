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
!SUBROUTINE closest_cable_segment
!
! NAME
!     SUBROUTINE closest_cable_segment
!
! DESCRIPTION
!
!     Search through the defined cable segments and return the closest segment which
!     contains the specified cable    
!
! COMMENTS
!     
!  
!
! HISTORY
!
!     started 21/09/2012 CJS
!
!
SUBROUTINE closest_cable_segment(point_number,cable_number,segment_number,segment_point)

USE TLM_general
USE geometry_types
USE geometry
USE Cables
USE constants

IMPLICIT NONE

  integer		:: point_number
  integer		:: cable_number
  integer		:: segment_number
  type(cell_point)	:: segment_point

! local variables

  integer	:: segment
  integer	:: cable
  integer	:: n_cables_found
  type(xyz)	:: point_coordinate
  type(xyz)	:: segment_coordinate1
  type(xyz)	:: segment_coordinate2
  real*8	:: dist1,dist2
  real*8	:: min_dist
  integer	:: min_segment
  type(cell_point)	:: min_segment_point
  
  logical	:: output_cable_found

! function variables 
 
  real*8	::xyz_distance

! START

  CALL write_line('CALLED: closest_cable_segment',0,output_to_screen_flag)

  output_cable_found=.FALSE.
  min_dist=1D30
  
! get the xyz coordinates of the required output point from the cell%i,cell%j,cell%k
  CALL get_cell_centre_coordinate(problem_points(point_number)%cell,point_coordinate)
    
  do segment=1,n_bundle_segments

! check whether the required cable exists in this bundle segment
    n_cables_found=0
    do cable=1,bundle_segment_list(segment)%n_cables
      if (bundle_segment_list(segment)%cable_list(cable).eq.cable_number) n_cables_found=n_cables_found+1
    end do
    
    if (n_cables_found.ne.0) then

! get segment end point coordinates  
      CALL get_cell_point_coordinate(bundle_segment_list(segment)%cable_segment%segment_point(1),segment_coordinate1)
      CALL get_cell_point_coordinate(bundle_segment_list(segment)%cable_segment%segment_point(2),segment_coordinate2)
    
      dist1=xyz_distance(point_coordinate,segment_coordinate1)
      if (dist1.lt.min_dist) then
        min_dist=dist1
	min_segment=segment
	min_segment_point=bundle_segment_list(segment)%cable_segment%segment_point(1)
        output_cable_found=.TRUE.
      end if
      
      dist2=xyz_distance(point_coordinate,segment_coordinate2)
      if (dist2.lt.min_dist) then
        min_dist=dist2
	min_segment=segment
	min_segment_point=bundle_segment_list(segment)%cable_segment%segment_point(2)
        output_cable_found=.TRUE.
      end if

    end if ! the required cable exists in this bundle segment
    
  end do ! next bundle segment
  
  if (.NOT.output_cable_found) GOTO 9000
  
  segment_number=min_segment
  segment_point=min_segment_point

  CALL write_line('FINISHED: closest_cable_segment',0,output_to_screen_flag)
    
  RETURN
  
9000 CALL write_line('ERROR in closest_cable_segment',0,.TRUE.)
     CALL write_line('Segment not found',0,.TRUE.)
     CALL write_line_integer('Point number',point_number,0,.TRUE.)
     CALL write_line_integer('Cable number',cable_number,0,.TRUE.)
     STOP
  
END SUBROUTINE closest_cable_segment
