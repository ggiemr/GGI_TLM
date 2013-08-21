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
!SUBROUTINE build_point_geometry
!
! NAME
!     SUBROUTINE build_point_geometry
!
! DESCRIPTION
!     build_point_geometry:
!     The only action required here is to apply the transformation to point coordinates. 
!     points are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE build_point_geometry

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: point_number

! START

  CALL write_line('CALLED: Build_point_geometry',0,output_to_screen_flag)

  do point_number=1,n_points

! apply the transformation to each of the point points
    CALL apply_transformation(problem_points(point_number)%point,problem_points(point_number)%trans)

  end do ! next point number
  
  CALL plot_points()

  CALL write_line('FINISHED: Build_point_geometry',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE build_point_geometry
