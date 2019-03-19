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
!     SUBROUTINE write_GGI_TLM_input_file
!
! DESCRIPTION
!	
!     
! COMMENTS
!     
!
!
!
! HISTORY
!
!     started 14/3/19 CJS
!     
!

SUBROUTINE write_GGI_TLM_input_file

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i,ii
character(LEN=256) :: line

! START

  write(*,*)'CALLED: write_GGI_TLM_input_file'

! mesh dimensions

  write(20,'(A)')'Mesh_outer_boundary_dimension'
  write(20,'(6ES16.6)')xmin,xmax,ymin,ymax,zmin,zmax 

! cell size

  write(20,*)
  write(20,'(A)')'Mesh_cell_dimension'
  write(20,'(ES16.6)')dl

! absorbing boundary condition

  write(20,*)
  write(20,'(A)')'Outer_boundary_reflection_coefficient'
  write(20,'(A)')'0.0 0.0 0.0 0.0 0.0 0.0'

! surface geometry

  write(20,*)
  write(20,'(A)')'Surface_list'
  write(20,'(I6,A22)')n_surfaces,'  # Number of surfaces'
  
  do i=1,n_surfaces
    
    write(20,'(I6,A18)')i,'  # SURFACE NUMBER'
  
    if (surface_type(i).EQ.surface_type_stl) then
    
      write(20,'(A)')'stl_triangulated_surface'
      write(20,'(A)')trim(surface_filename(i))
      write(20,'(A)')'1.0 0.0'
      write(20,'(3ES16.6)')0.0,0.0,0.0
      write(20,'(3ES16.6)')0.0,0.0,surface_z_offset(i)
      
    else if (surface_type(i).EQ.surface_type_zplane) then
    
      write(20,'(A)')'zplane'
      write(20,'(6ES16.6)')(surface_parameters(i,ii),ii=1,6)
      write(20,'(3ES16.6)')0.0,0.0,0.0
      write(20,'(3ES16.6)')0.0,0.0,0.0
   
    else
    
      write(*,*)'ERROR: Unknown surface geometry type',surface_type_stl
      STOP 1
      
    end if
      
  end do ! next surface

! surface materials* START of GGI_TLM

  write(20,*)
  write(20,'(A)')'Surface_material_list'
  write(20,'(I6,A31)')n_surface_materials,'  # Number of surface materials'
  
  do i=1,n_surface_materials
  
    write(20,'(I6,A27)')i,'  # SURFACE MATERIAL NUMBER'
    
    if ( surface_material_type(i).EQ.surface_material_type_PEC ) then
      write(20,'(A3)')'PEC'
      
      if (i.EQ.1) then                     ! special case for conducting surfaces from gerber files
      
        write(20,'(I6,A22)')n_PEC_surfaces,'  # Number of surfaces'      
        do ii=1,n_PEC_surfaces
          write(20,'(I5)',ADVANCE='NO')PEC_surface_list(ii)
        end do
        write(20,'(A)')' # Surface list'      
        do ii=1,n_PEC_surfaces
          write(20,'(I2)',ADVANCE='NO')PEC_surface_orientation_list(ii)
        end do
        write(20,'(A)')' # Surface orientation list'
        
      else
       
        write(20,'(A)')'     1  # Number of surfaces'
        write(20,'(I6,A18)')surface_material_to_surface_list(i),'  # Surface number'
        write(20,'(A)')'     1  # Surface orientation list'
      
      end if
      
    else if ( surface_material_type(i).EQ.surface_material_type_SPICE ) then
    
      write(20,'(A5)')'SPICE'
      write(20,'(I6,A23)')SPICE_node_list(i),'  # Ngspice node number'
      write(20,'(A2,A18)')SPICE_port_direction_list(i),'  # Port direction'
      
      write(20,'(A)')'     1  # Number of surfaces'
      write(20,'(I6,A18)')surface_material_to_surface_list(i),'  # Surface number'
      write(20,'(A)')'     1  # Surface orientation list'
    
    else
    
      write(*,*)'ERROR: unknown surface material type:',surface_material_type(i)
      STOP 1
       
    end if
    
  end do ! next surface material

! volume geometry

! volume materials

! Additional text from the GGI_TLM_create_PCB_simulation_model input file

10  read(10,'(A)',END=1000,ERR=1000)line
  
    if( index(line,'* START of GGI_TLM').EQ.0 ) GOTO 10    

20  read(10,'(A)',END=1000,ERR=1000)line
  
    if ( index(line,'* END of GGI_TLM').EQ.0 ) then
    
      write(20,'(A)')trim(line)
      GOTO 20
    
    end if
    
1000 CONTINUE

  write(*,*)'FINISHED: write_GGI_TLM_input_file'


RETURN  
  
END SUBROUTINE write_GGI_TLM_input_file
