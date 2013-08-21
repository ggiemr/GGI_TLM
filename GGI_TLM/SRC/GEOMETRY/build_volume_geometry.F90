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
!SUBROUTINE build_volume_geometry
!
! NAME
!     SUBROUTINE build_volume_geometry
!
! DESCRIPTION
!     build_volume_geometry:
!
!     Create volume geometric entities from the defined volume type
!     The geometric entities consist of trinagulated volumes
!     Triangulated volumes are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE build_volume_geometry

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number

! START

  CALL write_line('CALLED: Build_volume_geometry',0,output_to_screen_flag)

  do volume_number=1,n_volumes

    if (problem_volumes(volume_number)%volume_type.EQ.volume_type_sphere) then
      
      CALL build_volume_sphere(volume_number)
      
    else if (problem_volumes(volume_number)%volume_type.EQ.volume_type_rectangular_block) then
      
      CALL build_volume_rectangular_block(volume_number)

    else if (problem_volumes(volume_number)%volume_type.EQ.volume_type_rectangular_block2) then
      
      CALL build_volume_rectangular_block2(volume_number)

    else if (problem_volumes(volume_number)%volume_type.EQ.volume_type_cylinder) then
      
      CALL build_volume_cylinder(volume_number)

    else if (problem_volumes(volume_number)%volume_type.EQ.volume_type_tet) then
      
      CALL build_volume_tet(volume_number)

    else if (problem_volumes(volume_number)%volume_type.EQ.volume_type_pyramid) then
      
      CALL build_volume_pyramid(volume_number)

    else if (problem_volumes(volume_number)%volume_type.EQ.volume_type_pyramid_ram) then
      
      CALL build_volume_pyramid_ram(volume_number)

    else if (problem_volumes(volume_number)%volume_type.EQ.volume_type_tet_mesh) then
      
      CALL build_volume_tet_mesh(volume_number)
       
    else ! volume type not yet defined
    
      GOTO 9000
      
    end if ! volume_type

  end do ! next volume number
  
  CALL plot_tet_volumes()

  CALL write_line('FINISHED: Build_volume_geometry',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error in build_volume_geometry:',0,.TRUE.)
     CALL write_line('volume type number not defined',0,.TRUE.)
     CALL write_line_integer('volume type number',problem_volumes(volume_number)%volume_type,0,.TRUE.)
     STOP

  
END SUBROUTINE build_volume_geometry
