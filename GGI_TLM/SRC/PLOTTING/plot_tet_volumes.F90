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
!SUBROUTINE plot_tet_volumes
!
! NAME
!     SUBROUTINE plot_tet_volumes
!
! DESCRIPTION
!     plot_tet_volumes:
!
!     All the tet volumes are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE plot_tet_volumes()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number

integer	:: number_of_tets

! START

  CALL write_line('CALLED: plot_tet_volumes',0,output_to_screen_flag)

  do volume_number=1,n_volumes

    number_of_tets=problem_volumes(volume_number)%number_of_tets

    if (number_of_tets.gt.0) then

! open and write triangulated volume to vtk format file
      CALL open_vtk_file(tet_volume_file_unit,tet_volume_file_extension,volume_number) 
      
      CALL write_volume_list_vtk(tet_volume_file_unit,	&
                                  number_of_tets,problem_volumes(volume_number)%tet_list)
      
      CALL close_vtk_file(tet_volume_file_unit) 
    
    end if ! number_of_tets.gt.0

  end do ! next volume number

  CALL write_line('FINISHED: plot_tet_volumes',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_tet_volumes
