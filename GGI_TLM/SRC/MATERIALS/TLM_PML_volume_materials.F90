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
! SUBROUTINE set_pml_volume_material_mesh
! SUBROUTINE calculate_pml_volume_material_filter_coefficients
! SUBROUTINE allocate_pml_volume_material_filter_data
! SUBROUTINE pml_volume_material_update
!
! NAME
!     set_pml_volume_material_mesh
!
! DESCRIPTION
!     
!     Count the number of cells of each material and construct the 
!     appropriate lists for the main solver
!     
! COMMENTS
!     Allocate a local mesh array and fill with PML material data. This ensures that
!     a cell doesn't get given more than one material property.
!     In overlapping volumes the last material assignment takes precidence
!
! HISTORY
!
!     started 5/09/2012 CJS
!     frequency warping included 12/02/2013 CJS
!
!
SUBROUTINE set_pml_volume_material_mesh

USE TLM_general
USE geometry
USE PML_module
USE mesh
USE cell_parameters
USE file_information

IMPLICIT NONE

! local variables

  integer	:: volume_number
  integer	:: cell,number_of_cells
  integer	:: total_number_of_PML_cells
  integer	:: total_number_of_PML_cells_all_procsesses
  integer	:: n_data
  integer	:: cx,cy,cz
  
  integer	:: i
  integer	:: cell_number

! START
  
  CALL write_line('CALLED: set_pml_volume_material_mesh',0,output_to_screen_flag)

! INITIALISE VOLUME MATERIALS  
    
  CALL write_line_integer('n_pml_volume_materials=',n_pml_volume_materials,0,output_to_screen_flag)

  if (rank.eq.0) then
  
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'#START OF PML DESCRIPTION'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Number of PML volumes=',n_pml_volumes
    write(info_file_unit,*)''
  
  end if
  
! loop over the PML surfaces and flag the distance into the PML in x,y and z directions
! also check for intersection with other volume materials

  total_number_of_PMLl_cells=0

  do volume_number=1,n_pml_volumes
  
    pml_face=pml_volume_to_face(n_pml_volumes)

    number_of_cells=pml_volumes(volume_number)%number_of_cells

    CALL write_line_integer('PML volume     =',volume_number  ,0,output_to_screen_flag)
    CALL write_line_integer('Number of cells=',number_of_cells,0,output_to_screen_flag)

    if (number_of_cells.gt.0) then
      
      total_number_of_PMLl_cells=total_number_of_material_cells+number_of_cells
      
      do cell=1,number_of_cells
    
	cx=problem_volumes(volume_number)%cell_list(cell)%cell%i
	cy=problem_volumes(volume_number)%cell_list(cell)%cell%j
	cz=problem_volumes(volume_number)%cell_list(cell)%cell%k
	
	if (rank.eq.cell_rank(cz)) then
	  
          if (pml_face.EQ.pml_face_xmin) then
          
            local_cell_PML(cx,cy,cz,1)=dist_xmin
            
          else if (pml_face.EQ.pml_face_xmax) then
          
            local_cell_PML(cx,cy,cz,1)=dist_xmax
            
          else if (pml_face.EQ.pml_face_ymin) then
          
            local_cell_PML(cx,cy,cz,2)=dist_ymin
            
          else if (pml_face.EQ.pml_face_ymax) then
          
            local_cell_PML(cx,cy,cz,2)=dist_ymax
            
          else if (pml_face.EQ.pml_face_zmin) then
          
            local_cell_PML(cx,cy,cz,3)=dist_zmin
            
          else if (pml_face.EQ.pml_face_zmax) then
          
            local_cell_PML(cx,cy,cz,3)=dist_zmax
            
	  end if
          
	end if
	
      end do ! next cell
    
    end if ! number_of_cells.gt.0

  end do ! next volume number

****** WORKING HERE *********
  
!  local_cell_PML(1:nx,1:ny,nzmin:nzmax,1:3)=0

  do material_number=1,n_pml_volume_materials
  
    total_number_of_material_cells=0
  
    if (rank.eq.0) then 
      write(info_file_unit,*)''
      write(info_file_unit,*)'volume material_number : ',material_number  

      if (pml_volume_material_list(material_number)%type.EQ.pml_volume_material_type_PEC) then    
        write(info_file_unit,*)'Material type: PEC'  
      else if (pml_volume_material_list(material_number)%type.EQ.pml_volume_material_type_PMC) then 
        write(info_file_unit,*)'Material type: PMC'
      else if (pml_volume_material_list(material_number)%type.EQ.pml_volume_material_type_DISPERSIVE) then
        write(info_file_unit,*)'Material type: DISPERSIVE'
        write(info_file_unit,*)'Material name: ',trim(pml_volume_material_list(material_number)%name)
      end if
      
      CALL write_line_integer('Number of geometric volumes=',	&
                             pml_volume_material_list(material_number)%n_volumes,info_file_unit,output_to_screen_flag)
			     
    end if
    
    CALL write_line_integer('volume material number=',material_number,0,output_to_screen_flag)
    CALL write_line_integer('Number of geometric volumes=',	&
                             pml_volume_material_list(material_number)%n_volumes,0,output_to_screen_flag)

! loop over the geometric volumes with this material type and initialise the TLM cells
    do i=1,pml_volume_material_list(material_number)%n_volumes
    
      volume_number=pml_volume_material_list(material_number)%volume_list(i)

      CALL write_line_integer('Geometric volume material number=',volume_number,0,output_to_screen_flag)
      
      problem_volumes(volume_number)%pml_volume_material_number=material_number
      
      number_of_cells=problem_volumes(volume_number)%number_of_cells

      CALL write_line_integer('Number of cells=',number_of_cells,0,output_to_screen_flag)
      
      total_number_of_material_cells=total_number_of_material_cells+number_of_cells
      
      do cell=1,number_of_cells
    
	cx=problem_volumes(volume_number)%cell_list(cell)%cell%i
	cy=problem_volumes(volume_number)%cell_list(cell)%cell%j
	cz=problem_volumes(volume_number)%cell_list(cell)%cell%k
	
	if (rank.eq.cell_rank(cz)) then
	
          local_cell_material(cx,cy,cz)=material_number
	  
	end if
	
      end do ! next cell
    
    end do ! next volume number

#if defined(MPI)

    n_data=1
    call MPI_REDUCE(total_number_of_material_cells,total_number_of_material_cells_all_procsesses , &
                    n_data,MPI_INTEGER, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
		    
#elif defined(SEQ)

    total_number_of_material_cells_all_procsesses=total_number_of_material_cells
    
#endif			

    if (rank.eq.0) then
      CALL write_line_integer('Total number of material cells=',	&
                              total_number_of_material_cells_all_procsesses,0,output_to_screen_flag)
      write(info_file_unit,'(A,I14)')'Total number of material cells=',total_number_of_material_cells_all_procsesses
    end if

  end do ! next volume material number
  
  if (rank.eq.0) then
    write(info_file_unit,*)''
    write(info_file_unit,*)'#END OF VOLUME MATERIAL DESCRIPTION'
  end if
  
  CALL write_line('FINISHED: set_pml_volume_material_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_pml_volume_material_mesh
