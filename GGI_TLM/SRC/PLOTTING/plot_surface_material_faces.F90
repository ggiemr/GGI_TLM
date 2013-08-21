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
!SUBROUTINE plot_surface_material_faces
!
! NAME
!     SUBROUTINE plot_surface_material_faces
!
! DESCRIPTION
!     plot_surface_material_faces:
!
!     All the triangulated surfaces are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/11/2012 CJS
!     revised 19/2/2013 CJS -plot material faces as included in the TLM simulation
!                            Some may get overwritten by other materials
!
!
SUBROUTINE plot_surface_material_faces()

USE TLM_general
USE mesh
USE geometry_types
USE geometry
USE TLM_surface_materials
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer :: surface_material_number
integer	:: number_of_surfaces,surface_number
integer	:: total_number_of_surfaces
integer	:: total_number_of_faces

type(cell_point),allocatable	:: local_face_list(:)

integer :: i,count,face
integer :: cx,cy,cz

! START

  CALL write_line('CALLED: plot_surface_material_faces',0,output_to_screen_flag)

! read the surface material number to view  
  write(*,*)
  write(*,*)'Enter the surface material number to view'
  read(*,*)surface_material_number
  
  if ( (surface_material_number.le.0).OR.(surface_material_number.gt.n_surface_materials) ) then
    write(*,*)'surface_material_number is outside the available range'
    write(*,*)'Number of surface materials=',n_surface_materials
    RETURN
  end if
  
! work out the number of faces to write by looping over the surface list
  
  number_of_surfaces=surface_material_list(surface_material_number)%n_surfaces
  
  total_number_of_faces=0
  
  do cx=1,nx
    do cy=1,ny
      do cz=1,nz
        do face=1,3
	  if (abs(local_surface_material(cx,cy,cz,face)).EQ.surface_material_number) then
	    total_number_of_faces=total_number_of_faces+1
	  end if
	end do
      end do
    end do
  end do
  
    
  if (total_number_of_faces.EQ.0) then
    write(*,*)'Total number of faces of this material is 0'
    RETURN
  else
    write(*,*)'Total number of faces of this material is :',total_number_of_faces
  end if
  
! Allocate data to construct a local_face_list

  ALLOCATE( local_face_list(1:total_number_of_faces ) )
  
! fill the local_face_list  
  count=0
  do cx=1,nx
    do cy=1,ny
      do cz=1,nz
        do face=1,3
	
	  if (abs(local_surface_material(cx,cy,cz,face)).EQ.surface_material_number) then
            count=count+1
            local_face_list(count)%cell%i=cx
            local_face_list(count)%cell%j=cy
            local_face_list(count)%cell%k=cz
            local_face_list(count)%point=face
	  end if
	  
	end do
      end do
    end do
  end do
      
! open and write surface mesh to vtk format file
  CALL open_vtk_file(surface_material_faces_file_unit,surface_material_faces_file_extension,surface_material_number) 
      
  CALL write_surface_mesh_list_vtk(surface_material_faces_file_unit,	&
                                  total_number_of_faces,local_face_list)
      
  CALL close_vtk_file(surface_material_faces_file_unit) 

  
  DEALLOCATE( local_face_list )

  CALL write_line('FINISHED: plot_surface_material_faces',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_surface_material_faces
!SUBROUTINE plot_surface_material_faces_old
!
! NAME
!     SUBROUTINE plot_surface_material_faces
!
! DESCRIPTION
!     plot_surface_material_faces:
!
!     All the triangulated surfaces are written to .vtk files for visualisation
!     
! COMMENTS
!     This plots materials as specified on a surface basis - it does not plot
!     the computational surface mesh
!
! HISTORY
!
!     started 5/11/2012 CJS
!
!
SUBROUTINE plot_surface_material_faces_old()

USE TLM_general
USE geometry_types
USE geometry
USE TLM_surface_materials
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer :: surface_material_number
integer	:: number_of_surfaces,surface_number
integer	:: total_number_of_surfaces
integer	:: total_number_of_faces

type(cell_point),allocatable	:: local_face_list(:)

integer :: i,count,face

! START

  CALL write_line('CALLED: plot_surface_material_faces',0,output_to_screen_flag)

! read the surface material number to view  
  write(*,*)
  write(*,*)'Enter the surface material number to view'
  read(*,*)surface_material_number
  
  if ( (surface_material_number.le.0).OR.(surface_material_number.gt.n_surface_materials) ) then
    write(*,*)'surface_material_number is outside the available range'
    write(*,*)'Number of surface materials=',n_surface_materials
    RETURN
  end if
  
! work out the number of faces to write by looping over the surface list
  
  number_of_surfaces=surface_material_list(surface_material_number)%n_surfaces
  
  total_number_of_faces=0
  
  do i=1,number_of_surfaces
    surface_number=surface_material_list(surface_material_number)%surface_list(i)
    total_number_of_faces=total_number_of_faces+problem_surfaces(surface_number)%number_of_faces
  end do
  
  if (total_number_of_faces.EQ.0) then
    write(*,*)'Total number of faces of this material is 0'
    RETURN
  end if
  
! Allocate data to construct a local_face_list

  ALLOCATE( local_face_list(1:total_number_of_faces ) )
  
! fill the local_face_list  
  count=0
  do i=1,number_of_surfaces
    surface_number=surface_material_list(surface_material_number)%surface_list(i)
    do face=1,problem_surfaces(surface_number)%number_of_faces
      count=count+1
      local_face_list(count)=problem_surfaces(surface_number)%face_list(face)
    end do
  end do
      
! open and write surface mesh to vtk format file
  CALL open_vtk_file(surface_material_faces_file_unit,surface_material_faces_file_extension,surface_material_number) 
      
  CALL write_surface_mesh_list_vtk(surface_material_faces_file_unit,	&
                                  total_number_of_faces,local_face_list)
      
  CALL close_vtk_file(surface_material_faces_file_unit) 

  
  DEALLOCATE( local_face_list )

  CALL write_line('FINISHED: plot_surface_material_faces',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_surface_material_faces_old
