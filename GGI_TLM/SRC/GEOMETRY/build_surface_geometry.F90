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
!SUBROUTINE build_surface_geometry
!
! NAME
!     SUBROUTINE build_surface_geometry
!
! DESCRIPTION
!     build_surface_geometry:
!
!     Create surface geometric entities from the defined surface type
!     The geometric entities consist of triangulated surfaces
!     Triangulated surfaces are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE build_surface_geometry

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: surface_number

! START

  CALL write_line('CALLED: Build_surface_geometry',0,output_to_screen_flag)

  do surface_number=1,n_surfaces

    if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_sphere) then
      
      CALL build_surface_sphere(surface_number)
      
    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_rectangular_block) then
      
      CALL build_surface_rectangular_block(surface_number)

    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_rectangular_block2) then
      
      CALL build_surface_rectangular_block2(surface_number)

    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_cylinder) then
      
      CALL build_surface_cylinder(surface_number)

    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_rectangle) then
      
      CALL build_surface_rectangle(surface_number)
      
    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_xplane) then
      
      CALL build_surface_xplane(surface_number)
      
    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_yplane) then
      
      CALL build_surface_yplane(surface_number)
      
    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_zplane) then
      
      CALL build_surface_zplane(surface_number)

    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_circle) then
      
      CALL build_surface_circle(surface_number)

    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_quad) then
      
      CALL build_surface_quad(surface_number)

    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_triangle) then
      
      CALL build_surface_triangle(surface_number)
      
    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_triangulated_surface) then
      
      CALL build_surface_triangulated_surface(surface_number)
      
    else if (problem_surfaces(surface_number)%surface_type.EQ.surface_type_vtk_triangulated_surface) then
      
      CALL build_surface_vtk_triangulated_surface(surface_number)
       
    else ! surface type not yet defined
    
      GOTO 9000
      
    end if ! surface_type

  end do ! next surface number
  
  CALL plot_triangulated_surfaces()

  CALL write_line('FINISHED: Build_surface_geometry',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error in build_surface_geometry:',0,.TRUE.)
     CALL write_line('Surface type number not defined',0,.TRUE.)
     CALL write_line_integer('Surface type number',problem_surfaces(surface_number)%surface_type,0,.TRUE.)
     STOP

  
END SUBROUTINE build_surface_geometry
