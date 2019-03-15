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
!     SUBROUTINE add_components
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

SUBROUTINE add_components

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i,ii

real*8  :: sx,sy,sz
integer :: ix,iy,iz

integer :: component_number
character(LEN=256) :: line

! START

  read(10,*)n_lumped_components
  
  write(*,*)'Number of lumped components:',n_lumped_components
  
  do i=1,n_lumped_components
  
    read(10,*)component_number
  
    write(*,*)'Reading component number:',component_number
    
    if (component_number.NE.i) then
      write(*,*)'ERROR: components must be numbered in order, component number',i
      STOP 1
    end if
  
    read(10,'(A)')line  ! read the component type
    
    if (index(line,'one_port_model').NE.0) then
      lumped_component_type(i)=component_type_one_port_model
    else
      write(*,*)'ERROR: unknown component type:',trim(line)
      STOP 1
    end if
    
    read(10,*)ngspice_n_ports(i)
    
    write(*,*)'Number of ports',ngspice_n_ports(i)
    
    ngspice_n_nodes(i)=ngspice_n_ports(i)        ! number of ngspice nodes equal to number of ngspice ports
    
    write(*,*)'Number of ngspice nodes:',ngspice_n_nodes(i)
    
    tot_n_ngspice_nodes=tot_n_ngspice_nodes+ngspice_n_nodes(i)
    
    do ii=1,ngspice_n_nodes(i)
    
      read(10,*)ngspice_node_list(i,ii)
    
      write(*,*)'node',ii,' ngspice node number:',ngspice_node_list(i,ii)
    
    end do
    
    ngspice_n_terminals(i)=ngspice_n_ports(i)+1  ! number of terminals equal to number of ngspice ports +1
    
    write(*,*)'Number of terminals:',ngspice_n_terminals(i)

! read the coordinates of the terminals, also calculate the average position to determine where the active element will be
    sx=0d0
    sy=0d0
    sz=0d0
    
    do ii=1,ngspice_n_terminals(i)
    
      read(10,*)ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),ngspice_terminal_list(i,ii,3)
      
      write(*,*)'terminal',ii,' Coords:',ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),ngspice_terminal_list(i,ii,3)
      
      sx=sx+ngspice_terminal_list(i,ii,1)
      sy=sy+ngspice_terminal_list(i,ii,2)
      sz=sz+ngspice_terminal_list(i,ii,3)
    
    end do
    
    if (ngspice_n_terminals(i).NE.0) then
    
      sx=sx/dble(ngspice_n_terminals(i))
      sy=sy/dble(ngspice_n_terminals(i))
      sz=sz/dble(ngspice_n_terminals(i))
    
    else
      write(*,*)'ERROR calculating average position of component terminals'
      STOP 1
    end if
    
    read(10,*)z_position(i)
    
! Add the active element, single TLM face patch to the surface list
  
    n_surfaces=n_surfaces+1
    surface_type(n_surfaces)=surface_type_zplane
    
! calculate the position of the centre of the active element face
    sz=z_position(i)
    
! calculate the coordinates of the centre of the cell where the connecting wires will end

    CALL get_TLM_cell_from_coordinate(sx,sy,sz,ix,iy,iz)
    CALL get_TLM_cell_centre_coordinate(ix,iy,iz,sx,sy,sz)
    
    surface_parameters(n_surfaces,1)=sx-dl/2d0
    surface_parameters(n_surfaces,2)=sy-dl/2d0
    surface_parameters(n_surfaces,3)=sz-dl/2d0  ! zmin face
    surface_parameters(n_surfaces,4)=sx+dl/2d0
    surface_parameters(n_surfaces,5)=sy+dl/2d0
    surface_parameters(n_surfaces,6)=sz-dl/2d0  ! zmin face

! Add the active element material of type 'SPICE' to the surface material list

    n_surface_materials=n_surface_materials+1  
    
    surface_material_type(n_surface_materials)=surface_material_type_SPICE  
    
    surface_material_to_surface_list(n_surface_materials)=n_surfaces
    SPICE_node_list(n_surface_materials)=ngspice_node_list(i,1)      ! ****** ASSUME ONLY SINGLE NGSPICE NODE FOR NOW ******
    SPICE_port_direction_list(n_surface_materials)='-y'

! Set the cells from each terminal to the active element cell to PEC
    do ii=1,ngspice_n_terminals(i)
    
      write(*,*)'Add terminal connection cells:'
      write(*,'(3es16.6)')ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),ngspice_terminal_list(i,ii,3)
      write(*,'(3es16.6)') sx,sy,sz
    
      CALL set_terminal_connection_cells(ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),  &
                                         ngspice_terminal_list(i,ii,3),sx,sy,sz)
          
    end do

! Set the active element cell to free space - the active element is set elsewhere
    CALL get_TLM_cell_from_coordinate(sx,sy,sz,ix,iy,iz)
    material_mesh(centre,ix,iy,iz)=0
  
  end do ! next component
  
! Add a surface for all the PEC electrical connections between the active elements and their terminals
  
  n_surfaces=n_surfaces+1
  surface_type(n_surfaces)=surface_type_stl
  surface_filename(n_surfaces)=terminal_connection_geometry_filename
  surface_z_offset(n_surfaces)=0d0

! Add the active element material of type 'PEC' to the surface material list

  n_surface_materials=n_surface_materials+1  
  surface_material_type(n_surface_materials)=surface_material_type_PEC  
  surface_material_to_surface_list(n_surface_materials)=n_surfaces
  

RETURN  
  
END SUBROUTINE add_components
