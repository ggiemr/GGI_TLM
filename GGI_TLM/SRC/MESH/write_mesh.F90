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
!SUBROUTINE write_mesh
!
! NAME
!     SUBROUTINE write_mesh
!
! DESCRIPTION
!     write_mesh:
!
!     write volume, surface, line and point meshes to file
!
!     Include mesh copying required by the periodic boundary condition implementation 
!     using the method of Lee and Smith, 
!     "An alternative approach for implementing periodic boundary conditions in the FDTD method using multiple unit cells,"
!     IEEE trans AP vol 54, no2, 2006 pp 698-705
! 
!     
! COMMENTS
!     We only copy the surface and volume meshes for the periodic boundary stuff as wires are not yet included 
!     and output (i.e. on points) is only required in one of the four cells
!
! HISTORY
!
!     started 12/08/2012 CJS
!     periodic boundary mesh copy 5/03/2014 CJS
!             1/10/2019 CJS  Add PML
!
!
SUBROUTINE write_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE PML_module
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number
integer	:: cell,number_of_cells

integer	:: surface_number
integer	:: face,number_of_faces

integer	:: line_number
integer	:: segment,number_of_cell_segments

integer	:: point_number

integer :: nx2,ny2

integer :: i

! START

  CALL write_line('CALLED: write_mesh',0,output_to_screen_flag)

! Open mesh file

  CALL open_file(mesh_file_unit,mesh_file_extension)

! Write general mesh parameters

!  write(mesh_file_unit,*)periodic_boundary,' periodic_boundary'

  if (periodic_boundary) then
! double the mesh size in the x and y directions
    nx2=nx
    ny2=ny
    nx=nx*2
    ny=ny*2
    mesh_xmax=mesh_xmin+nx*dl
    mesh_ymax=mesh_ymin+ny*dl
  end if

  write(mesh_file_unit,*)dl,' dl'
  write(mesh_file_unit,*)nx,ny,nz,' nx ny nz'
  write(mesh_file_unit,*)mesh_xmin,mesh_xmax,' mesh_xmin,mesh_xmax'
  write(mesh_file_unit,*)mesh_ymin,mesh_ymax,' mesh_ymin,mesh_ymax'
  write(mesh_file_unit,*)mesh_zmin,mesh_zmax,' mesh_zmin,mesh_zmax'
  
  write(mesh_file_unit,*)n_volumes,' n_volumes'

  do volume_number=1,n_volumes

    number_of_cells=problem_volumes(volume_number)%number_of_cells
    
    if (.NOT.periodic_boundary) then
    
      write(mesh_file_unit,*)number_of_cells,' number_of_cells in volume=',volume_number
    
      do cell=1,number_of_cells
     
        write(mesh_file_unit,*)	problem_volumes(volume_number)%cell_list(cell)%cell%i,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%j,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%k,     &
				problem_volumes(volume_number)%cell_list(cell)%point
				
      end do ! next cell
     
    else ! each volume cell gets copied so we have 4 cells in the final model
    
      write(mesh_file_unit,*)number_of_cells*4,' number_of_cells in volume=',volume_number

! original cell    
      do cell=1,number_of_cells     
        write(mesh_file_unit,*)	problem_volumes(volume_number)%cell_list(cell)%cell%i,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%j,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%k,     &
				problem_volumes(volume_number)%cell_list(cell)%point				
      end do ! next cell

! cell +nx2
      do cell=1,number_of_cells     
        write(mesh_file_unit,*)	problem_volumes(volume_number)%cell_list(cell)%cell%i+nx2,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%j,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%k,     &
				problem_volumes(volume_number)%cell_list(cell)%point				
      end do ! next cell

! cell +ny2
      do cell=1,number_of_cells     
        write(mesh_file_unit,*)	problem_volumes(volume_number)%cell_list(cell)%cell%i,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%j+ny2,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%k,     &
				problem_volumes(volume_number)%cell_list(cell)%point				
      end do ! next cell
   
! cell +nx2+ny2
      do cell=1,number_of_cells     
        write(mesh_file_unit,*)	problem_volumes(volume_number)%cell_list(cell)%cell%i+nx2,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%j+ny2,     &
				problem_volumes(volume_number)%cell_list(cell)%cell%k,     &
				problem_volumes(volume_number)%cell_list(cell)%point				
      end do ! next cell
   
    end if  ! periodic boundary
    
  end do ! next volume number
  
  write(mesh_file_unit,*)n_surfaces,' n_surfaces'

  do surface_number=1,n_surfaces

    number_of_faces=problem_surfaces(surface_number)%number_of_faces
    
    if (.NOT.periodic_boundary) then
      write(mesh_file_unit,*)number_of_faces,' number_of_faces in surface=',surface_number
    else
      write(mesh_file_unit,*)number_of_faces*4,' number_of_faces in surface=',surface_number
    end if
    
    if (periodic_boundary) then
! the extent of the surface mesh is extended by the additional unit cells added
      problem_surfaces(surface_number)%mesh_xmax=problem_surfaces(surface_number)%mesh_xmax+nx2*dl
      problem_surfaces(surface_number)%mesh_ymax=problem_surfaces(surface_number)%mesh_ymax+ny2*dl
      problem_surfaces(surface_number)%mesh_cell_xmax=problem_surfaces(surface_number)%mesh_cell_xmax+nx2
      problem_surfaces(surface_number)%mesh_cell_ymax=problem_surfaces(surface_number)%mesh_cell_ymax+ny2
    end if

    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_xmin,        &
    			   problem_surfaces(surface_number)%mesh_xmax,' mesh_xmin,mesh_xmax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_ymin,        &
    			   problem_surfaces(surface_number)%mesh_ymax,' mesh_ymin,mesh_ymax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_zmin,        &
    			   problem_surfaces(surface_number)%mesh_zmax,' mesh_zmin,mesh_zmax'

    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_xmin,   &
    			   problem_surfaces(surface_number)%mesh_cell_xmax,' mesh_cell_xmin,mesh_cell_xmax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_ymin,   &
    			   problem_surfaces(surface_number)%mesh_cell_ymax,' mesh_cell_ymin,mesh_cell_ymax'
    write(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_zmin,   &
                           problem_surfaces(surface_number)%mesh_cell_zmax,' mesh_cell_zmin,mesh_cell_zmax'
     
    if (.NOT.periodic_boundary) then
   
      do face=1,number_of_faces
        write(mesh_file_unit,*)	problem_surfaces(surface_number)%face_list(face)%cell%i,	&
				problem_surfaces(surface_number)%face_list(face)%cell%j,	&
				problem_surfaces(surface_number)%face_list(face)%cell%k,	&
				problem_surfaces(surface_number)%face_list(face)%point
      end do ! next face
    
    else ! periodic boundary
    
! original cell    
      do face=1,number_of_faces
        write(mesh_file_unit,*)	problem_surfaces(surface_number)%face_list(face)%cell%i,	&
				problem_surfaces(surface_number)%face_list(face)%cell%j,	&
				problem_surfaces(surface_number)%face_list(face)%cell%k,	&
				problem_surfaces(surface_number)%face_list(face)%point
      end do ! next face
    
! cell +nx2
      do face=1,number_of_faces
        write(mesh_file_unit,*)	problem_surfaces(surface_number)%face_list(face)%cell%i+nx2,	&
				problem_surfaces(surface_number)%face_list(face)%cell%j,	&
				problem_surfaces(surface_number)%face_list(face)%cell%k,	&
				problem_surfaces(surface_number)%face_list(face)%point
      end do ! next face
    
! cell +ny2
      do face=1,number_of_faces
        write(mesh_file_unit,*)	problem_surfaces(surface_number)%face_list(face)%cell%i,	&
				problem_surfaces(surface_number)%face_list(face)%cell%j+ny2,	&
				problem_surfaces(surface_number)%face_list(face)%cell%k,	&
				problem_surfaces(surface_number)%face_list(face)%point
      end do ! next face
    
! cell +nx2 +ny2
      do face=1,number_of_faces
        write(mesh_file_unit,*)	problem_surfaces(surface_number)%face_list(face)%cell%i+nx2,	&
				problem_surfaces(surface_number)%face_list(face)%cell%j+ny2,	&
				problem_surfaces(surface_number)%face_list(face)%cell%k,	&
				problem_surfaces(surface_number)%face_list(face)%point
      end do ! next face
    
    end if  ! periodic boundary
    
  end do ! next surface number
  
  write(mesh_file_unit,*)n_lines,' n_lines'

  do line_number=1,n_lines

    number_of_cell_segments=problem_lines(line_number)%number_of_cell_segments
    
    write(mesh_file_unit,*)number_of_cell_segments,' number_of_cell_segments in line=',line_number
    
    do segment=1,number_of_cell_segments
      write(mesh_file_unit,*)	problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%i,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%j,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%cell%k,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(1)%point,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%i,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%j,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%cell%k,	&
				problem_lines(line_number)%cell_segment_list(segment)%segment_point(2)%point
    end do ! next cell segment
    
  end do ! next line number
  
  write(mesh_file_unit,*)n_points,' n_points'

  do point_number=1,n_points

    write(mesh_file_unit,*)	problem_points(point_number)%point%x,	&
				problem_points(point_number)%point%y,	&
				problem_points(point_number)%point%z
    write(mesh_file_unit,*)	problem_points(point_number)%cell%i,	&
				problem_points(point_number)%cell%j,	&
				problem_points(point_number)%cell%k
    write(mesh_file_unit,*)	problem_points(point_number)%face%cell%i,	&
				problem_points(point_number)%face%cell%j,	&
				problem_points(point_number)%face%cell%k,	&
				problem_points(point_number)%face%point
    
  end do ! next point number
  
  if ( (n_pml_volumes.NE.0).AND.periodic_boundary ) then
    write(*,*)'ERROR: We cannot use the PML with periodic boundaries'
    STOP 1
  end if
  
  write(mesh_file_unit,*)n_pml_volumes,' n_pml_volumes'

  do volume_number=1,n_pml_volumes

    number_of_cells=pml_volumes(volume_number)%number_of_cells
        
      write(mesh_file_unit,*)number_of_cells,' number_of_cells in PML volume=',volume_number
    
      do cell=1,number_of_cells
     
        write(mesh_file_unit,*)	pml_volumes(volume_number)%cell_list(cell)%cell%i,     &
				pml_volumes(volume_number)%cell_list(cell)%cell%j,     &
				pml_volumes(volume_number)%cell_list(cell)%cell%k,     &
				pml_volumes(volume_number)%cell_list(cell)%point
				
      end do ! next cell
    
  end do ! next volume number
  
  write(mesh_file_unit,*)'pml_volume_to_face'
  write(mesh_file_unit,*)(pml_volume_to_face(i),i=1,6)
  
  write(mesh_file_unit,*)'pml thicknesses on xmin, xmax, ymin, ymax, zmin and zmax boundaries'
  write(mesh_file_unit,*)pml_txmin,pml_txmax
  write(mesh_file_unit,*)pml_tymin,pml_tymax
  write(mesh_file_unit,*)pml_tzmin,pml_tzmax
  
  write(mesh_file_unit,*)'pml_r'
  write(mesh_file_unit,*)pml_r
  
  write(mesh_file_unit,*)'pml_order'
  write(mesh_file_unit,*)pml_order
  
  CALL close_file(mesh_file_unit)

  CALL write_line('FINISHED: write_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE write_mesh
