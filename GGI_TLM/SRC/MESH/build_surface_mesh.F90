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
!SUBROUTINE build_surface_mesh
!
! NAME
!     SUBROUTINE build_surface_mesh
!
! DESCRIPTION
!     build_surface_mesh: Create a surface mesh array and put meshed triangles onto it
!     then count the cell surfaces and write into the problem_surfaces(surface_number)%face_list
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE build_surface_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: surface_number

integer,allocatable	:: surface_mesh(:,:,:,:)

integer			:: number_of_triangles

integer			:: triangle_number
integer			:: face_count

integer			:: ix,iy,iz,face

type(xyz)		:: xyz_point

! START

  CALL write_line('CALLED: build_surface_mesh',0,output_to_screen_flag)
  
  ALLOCATE( surface_mesh(1:nx,1:ny,1:nz,1:6) )

  do surface_number=1,n_surfaces
  
! reset surface mesh
    surface_mesh(1:nx,1:ny,1:nz,1:6)=0
    
    number_of_triangles=problem_surfaces(surface_number)%number_of_triangles

! Mesh each triangle
    do triangle_number=1,number_of_triangles    
    
      if (new_mesh_generation) then
      
        CALL mesh_triangle_new(surface_mesh,problem_surfaces(surface_number)%triangle_list(triangle_number))
      
      else
      
        CALL mesh_triangle(surface_mesh,problem_surfaces(surface_number)%triangle_list(triangle_number))
      
      end if
      
    end do ! next triangle
    
    
! count cell faces
    face_count=0

    do face=1,6
      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
	  
	    if (surface_mesh(ix,iy,iz,face).NE.0) face_count=face_count+1
	    
	  end do ! ix
        end do ! iy
      end do ! iz
    end do ! next face
    
!    write(*,*)'Number of faces=',face_count   
!    write(*,*)'FINISHED MESHING TRIANGLE'
!    STOP   

! Allocate memory for the mesh

    problem_surfaces(surface_number)%number_of_faces=face_count
    
    if (face_count.gt.0) then
    
      ALLOCATE( problem_surfaces(surface_number)%face_list(1:face_count) )

! loop over mesh surfaces and fill the face list, also work out the min and max dimensions of
! the mesh in each direction

      face_count=0
      
      problem_surfaces(surface_number)%mesh_xmin=mesh_xmax
      problem_surfaces(surface_number)%mesh_xmax=mesh_xmin
      problem_surfaces(surface_number)%mesh_ymin=mesh_ymax
      problem_surfaces(surface_number)%mesh_ymax=mesh_ymin
      problem_surfaces(surface_number)%mesh_zmin=mesh_zmax
      problem_surfaces(surface_number)%mesh_zmax=mesh_zmin
      
      problem_surfaces(surface_number)%mesh_cell_xmin=nx
      problem_surfaces(surface_number)%mesh_cell_xmax=0
      problem_surfaces(surface_number)%mesh_cell_ymin=ny
      problem_surfaces(surface_number)%mesh_cell_ymax=0
      problem_surfaces(surface_number)%mesh_cell_zmin=nz
      problem_surfaces(surface_number)%mesh_cell_zmax=0
 
      do face=1,6
        do iz=1,nz
          do iy=1,ny
            do ix=1,nx
	  
	      if (surface_mesh(ix,iy,iz,face).NE.0) then
	      
	        face_count=face_count+1
	        problem_surfaces(surface_number)%face_list(face_count)%cell%i=ix
	        problem_surfaces(surface_number)%face_list(face_count)%cell%j=iy
	        problem_surfaces(surface_number)%face_list(face_count)%cell%k=iz
	        problem_surfaces(surface_number)%face_list(face_count)%point=face
		
		CALL get_cell_point_coordinate(problem_surfaces(surface_number)%face_list(face_count),xyz_point)
		
                problem_surfaces(surface_number)%mesh_xmin=min(problem_surfaces(surface_number)%mesh_xmin,xyz_point%x)
                problem_surfaces(surface_number)%mesh_xmax=max(problem_surfaces(surface_number)%mesh_xmax,xyz_point%x)
                problem_surfaces(surface_number)%mesh_ymin=min(problem_surfaces(surface_number)%mesh_ymin,xyz_point%y)
                problem_surfaces(surface_number)%mesh_ymax=max(problem_surfaces(surface_number)%mesh_ymax,xyz_point%y)
                problem_surfaces(surface_number)%mesh_zmin=min(problem_surfaces(surface_number)%mesh_zmin,xyz_point%z)
                problem_surfaces(surface_number)%mesh_zmax=max(problem_surfaces(surface_number)%mesh_zmax,xyz_point%z)
		
                problem_surfaces(surface_number)%mesh_cell_xmin=min(problem_surfaces(surface_number)%mesh_cell_xmin,ix)
                problem_surfaces(surface_number)%mesh_cell_xmax=max(problem_surfaces(surface_number)%mesh_cell_xmax,ix)
                problem_surfaces(surface_number)%mesh_cell_ymin=min(problem_surfaces(surface_number)%mesh_cell_ymin,iy)
                problem_surfaces(surface_number)%mesh_cell_ymax=max(problem_surfaces(surface_number)%mesh_cell_ymax,iy)
                problem_surfaces(surface_number)%mesh_cell_zmin=min(problem_surfaces(surface_number)%mesh_cell_zmin,iz)
                problem_surfaces(surface_number)%mesh_cell_zmax=max(problem_surfaces(surface_number)%mesh_cell_zmax,iz)
		
              end if ! this is a surface face
	      
	    end do ! ix
          end do ! iy
        end do ! iz
      end do ! next face
      
    end if ! number of mesh faces .gt.0
    
  end do ! next surface number
  
  DEALLOCATE( surface_mesh )

  CALL write_line('FINISHED: build_surface_mesh',0,output_to_screen_flag)

  RETURN
  
  
END SUBROUTINE build_surface_mesh
