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
!
! NAME
!     SUBROUTINE write_volume_mesh_list_vtk
!
! DESCRIPTION
!     write_volume_mesh_list_vtk:
!
!     Write cell list to vtk format file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!     29/04/2013 Only write cell faces where there is a material discontinuity i.e. the external 
!                surface of the volume CJS
!
!
SUBROUTINE write_volume_mesh_list_vtk(file_unit,number_of_cells,cell_list)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE cell_parameters
USE constants

IMPLICIT NONE

integer	:: file_unit
integer	:: number_of_cells

type(cell_point)		:: cell_list(1:number_of_cells)

! local variables

integer 			:: cell_number

integer,allocatable		:: local_mesh(:,:,:)

integer				:: ix,iy,iz
integer 			:: loop

integer				:: number_of_faces
type(cell_point),allocatable	:: face_list(:)

! START

  CALL write_line('CALLED: Write_volume_mesh_list',0,output_to_screen_flag)



! Create a volume mesh and fill it with the volume material

  ALLOCATE ( local_mesh(0:nx+1,0:ny+1,0:nz+1) )
  
  local_mesh(0:nx+1,0:ny+1,0:nz+1) =0
    
  do cell_number=1,number_of_cells
  
    ix=cell_list(cell_number)%cell%i
    iy=cell_list(cell_number)%cell%j
    iz=cell_list(cell_number)%cell%k
    local_mesh(ix,iy,iz)=1
    
  end do! next cell 

! Find all the external surfaces of the volume material and fill the face list

  do loop=1,2
  
    number_of_faces=0

    do iz=1,nz
      do iy=1,ny	
	do ix=1,nx

! faces normal to x
          if ( (local_mesh(ix,iy,iz).EQ.1).AND.(local_mesh(ix-1,iy,iz).EQ.0) ) then ! include xmin face
            number_of_faces=number_of_faces+1
	    if (loop.eq.2) then
	      face_list(number_of_faces)%cell%i=ix
	      face_list(number_of_faces)%cell%j=iy
	      face_list(number_of_faces)%cell%k=iz
	      face_list(number_of_faces)%point=face_xmin
	    end if ! loop=2
          end if ! include xmin face

          if ( (local_mesh(ix,iy,iz).EQ.1).AND.(local_mesh(ix+1,iy,iz).EQ.0) ) then ! include xmax face
            number_of_faces=number_of_faces+1
	    if (loop.eq.2) then
	      face_list(number_of_faces)%cell%i=ix
	      face_list(number_of_faces)%cell%j=iy
	      face_list(number_of_faces)%cell%k=iz
	      face_list(number_of_faces)%point=face_xmax
	    end if ! loop=2
          end if ! include xmax face

! faces normal to y
          if ( (local_mesh(ix,iy,iz).EQ.1).AND.(local_mesh(ix,iy-1,iz).EQ.0) ) then ! include ymin face
            number_of_faces=number_of_faces+1
	    if (loop.eq.2) then
	      face_list(number_of_faces)%cell%i=ix
	      face_list(number_of_faces)%cell%j=iy
	      face_list(number_of_faces)%cell%k=iz
	      face_list(number_of_faces)%point=face_ymin
	    end if ! loop=2
          end if ! include ymin face

          if ( (local_mesh(ix,iy,iz).EQ.1).AND.(local_mesh(ix,iy+1,iz).EQ.0) ) then ! include ymax face
            number_of_faces=number_of_faces+1
	    if (loop.eq.2) then
	      face_list(number_of_faces)%cell%i=ix
	      face_list(number_of_faces)%cell%j=iy
	      face_list(number_of_faces)%cell%k=iz
	      face_list(number_of_faces)%point=face_ymax
	    end if ! loop=2
          end if ! include ymax face

! faces normal to z
          if ( (local_mesh(ix,iy,iz).EQ.1).AND.(local_mesh(ix,iy,iz-1).EQ.0) ) then ! include zmin face
            number_of_faces=number_of_faces+1
	    if (loop.eq.2) then
	      face_list(number_of_faces)%cell%i=ix
	      face_list(number_of_faces)%cell%j=iy
	      face_list(number_of_faces)%cell%k=iz
	      face_list(number_of_faces)%point=face_zmin
	    end if ! loop=2
          end if ! include zmin face

          if ( (local_mesh(ix,iy,iz).EQ.1).AND.(local_mesh(ix,iy,iz+1).EQ.0) ) then ! include zmax face
            number_of_faces=number_of_faces+1
	    if (loop.eq.2) then
	      face_list(number_of_faces)%cell%i=ix
	      face_list(number_of_faces)%cell%j=iy
	      face_list(number_of_faces)%cell%k=iz
	      face_list(number_of_faces)%point=face_zmax
	    end if ! loop=2
          end if ! include zmax face
	
	end do	! next ix
      end do  ! next iy
    end do  ! next iz
    
    if (loop.eq.1) then
      ALLOCATE ( face_list(1:number_of_faces) )
    end if
    
  end do

! plot the face list

  CALL write_surface_mesh_list_vtk(file_unit,number_of_faces,face_list)

  DEALLOCATE ( local_mesh )
  DEALLOCATE ( face_list )

  CALL write_line('FINISHED: Write_volume_mesh_list',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE write_volume_mesh_list_vtk
