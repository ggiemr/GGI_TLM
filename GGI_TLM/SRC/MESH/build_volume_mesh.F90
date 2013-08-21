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
!SUBROUTINE build_volume_mesh
!
! NAME
!     SUBROUTINE build_volume_mesh
!
! DESCRIPTION
!     build_volume_mesh:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
SUBROUTINE build_volume_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number

integer	:: number_of_tets
integer	:: tet_number
integer	:: cell_count

integer,allocatable	:: volume_mesh(:,:,:)

integer			:: ix,iy,iz

! START

  CALL write_line('CALLED: build_volume_mesh',0,output_to_screen_flag)
  
  ALLOCATE( volume_mesh(1:nx,1:ny,1:nz) )

  do volume_number=1,n_volumes

! reset volume mesh
    volume_mesh(1:nx,1:ny,1:nz)=0
    
    number_of_tets=problem_volumes(volume_number)%number_of_tets

    do tet_number=1,number_of_tets    
    
      if (new_mesh_generation) then
      
        CALL mesh_tet_new(volume_mesh,problem_volumes(volume_number)%tet_list(tet_number))
      
      else
      
        CALL mesh_tet(volume_mesh,problem_volumes(volume_number)%tet_list(tet_number))
       
      end if
     
    end do ! next tet
    
! count cells
    cell_count=0

    do iz=1,nz
      do iy=1,ny
    	do ix=1,nx
        
          if (volume_mesh(ix,iy,iz).NE.0) cell_count=cell_count+1
          
        end do ! ix
      end do ! iy
    end do ! iz

! Allocate memory for the mesh

    problem_volumes(volume_number)%number_of_cells=cell_count
    
    if (cell_count.gt.0) then
    
      ALLOCATE( problem_volumes(volume_number)%cell_list(1:cell_count) )

! loop over mesh allocating the cell list    
      cell_count=0
      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
	  
	    if (volume_mesh(ix,iy,iz).NE.0) then
	    
	      cell_count=cell_count+1
	      problem_volumes(volume_number)%cell_list(cell_count)%cell%i=ix
	      problem_volumes(volume_number)%cell_list(cell_count)%cell%j=iy
	      problem_volumes(volume_number)%cell_list(cell_count)%cell%k=iz
	      problem_volumes(volume_number)%cell_list(cell_count)%point=centre
	      
	    end if ! this cell is in the mesh
	    
	  end do ! ix
        end do ! iy
      end do ! iz

      
    end if ! number of mesh cells .gt.0
    
  end do ! next volume number
  
  DEALLOCATE( volume_mesh )

  CALL write_line('FINISHED: build_volume_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE build_volume_mesh
