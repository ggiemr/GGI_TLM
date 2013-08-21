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
!SUBROUTINE plot_triangulated_surfaces
!
! NAME
!     SUBROUTINE plot_triangulated_surfaces
!
! DESCRIPTION
!     plot_triangulated_surfaces:
!
!     All the triangulated surfaces are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE plot_triangulated_surfaces()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: surface_number

integer	:: number_of_triangles

! START

  CALL write_line('CALLED: Plot_triangulated_surfaces',0,output_to_screen_flag)

  do surface_number=1,n_surfaces

    number_of_triangles=problem_surfaces(surface_number)%number_of_triangles

    if (number_of_triangles.gt.0) then

! open and write triangulated surface to vtk format file
      CALL open_vtk_file(triangulated_surface_file_unit,triangulated_surface_file_extension,surface_number) 
      
      CALL write_surface_list_vtk(triangulated_surface_file_unit,	&
                                  number_of_triangles,problem_surfaces(surface_number)%triangle_list)
      
      CALL close_vtk_file(triangulated_surface_file_unit) 
    
    end if ! number_of_triangles.gt.0

  end do ! next surface number

  CALL write_line('FINISHED: Plot_triangulated_surfaces',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_triangulated_surfaces
