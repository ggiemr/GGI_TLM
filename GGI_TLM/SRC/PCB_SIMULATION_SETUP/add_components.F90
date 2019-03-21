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
!  I think that the port direction specification needs to be thought about...
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

integer :: i,ii,cell

real*8  :: sx,sy,sz        ! centre coordinates of device
integer :: ix,iy,iz        ! centre cell of device

real*8  :: sxmin,symin,szmin        ! minimum coordinates of device
real*8  :: sxmax,symax,szmax        ! maximum coordinates of device
real*8  :: lx,ly,lz        ! extent of device

integer :: dx,dy,dz        ! direction of multi-port device

real*8  :: sxp(5),syp(5),szp(5)     ! coordinates of two port device cell centres
integer :: ixp(5),iyp(5),izp(5)     ! centre cells of two port device cell centres

integer :: component_number
character(LEN=256) :: line

character(LEN=256) :: material_name

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
      
    else if (index(line,'two_port_model').NE.0) then
    
      lumped_component_type(i)=component_type_two_port_model
      
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
! also the x, y and z extent of the device

    sx=0d0
    sy=0d0
    sz=0d0
    
    sxmin=1d30
    symin=1d30
    szmin=1d30
    sxmax=-1d30
    symax=-1d30
    szmax=-1d30
    
    do ii=1,ngspice_n_terminals(i)
    
      read(10,*)ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),ngspice_terminal_list(i,ii,3)
      
      write(*,*)'terminal',ii,' Coords:',ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),ngspice_terminal_list(i,ii,3)
      
      sx=sx+ngspice_terminal_list(i,ii,1)
      sy=sy+ngspice_terminal_list(i,ii,2)
      sz=sz+ngspice_terminal_list(i,ii,3)
      
      sxmin=min(sxmin,ngspice_terminal_list(i,ii,1))
      symin=min(symin,ngspice_terminal_list(i,ii,2))
      szmin=min(szmin,ngspice_terminal_list(i,ii,3))
      sxmax=max(sxmax,ngspice_terminal_list(i,ii,1))
      symax=max(symax,ngspice_terminal_list(i,ii,2))
      szmax=max(szmax,ngspice_terminal_list(i,ii,3))
    
    end do
    
    if (ngspice_n_terminals(i).NE.0) then
    
      sx=sx/dble(ngspice_n_terminals(i))
      sy=sy/dble(ngspice_n_terminals(i))
      sz=sz/dble(ngspice_n_terminals(i))
    
    else
      write(*,*)'ERROR calculating average position of component terminals'
      STOP 1
    end if
    
    lx=sxmax-sxmin
    ly=symax-symin
    lz=szmax-szmin
    
    write(*,*)'min coordinate of component:'
    write(*,*)sxmin,symin,szmin
    write(*,*)'max coordinate of component:'
    write(*,*)sxmax,symax,szmax
    write(*,*)'extent of component:'
    write(*,*)lx,ly,lz
    
    read(10,*)z_position(i)
    
    if (ngspice_n_terminals(i).EQ.2) then
    
! **** TWO TERMINAL DEVICE ****
    
! calculate the position of the centre of the active element face
      sz=z_position(i)

! Add the active element, single TLM face patch to the surface list
  
      n_surfaces=n_surfaces+1
      surface_type(n_surfaces)=surface_type_zplane
    
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
      SPICE_node_list(n_surface_materials)=ngspice_node_list(i,1)     
      SPICE_port_direction_list(n_surface_materials)='-y'         ! ****** PORT DIRECTION TO BE GENERALISED ******

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
  
    else if (ngspice_n_terminals(i).EQ.3) then
    
! **** THREE TERMINAL DEVICE ****
    
! calculate the position of the centre of the device
      sz=z_position(i)
      
! work out the best direction for the internal layout of the device based on the largest dimension between terminals
      dz=0                 ! Assume the active elements are normal to z
      dx=0
      dy=0
      
      if (lx.GT.ly) then
        dx=1
      else
        dy=1
      end if 
      
      write(*,*)'Component direction dx,dy,dz:'
      write(*,*)dx,dy,dz
      
! Calculate the coordinates of the five cells constituting the terminals and active element cells
   
! calculate the coordinates of the centre cell where the central wire will end

      CALL get_TLM_cell_from_coordinate(sx,sy,sz,ix,iy,iz)
      CALL get_TLM_cell_centre_coordinate(ix,iy,iz,sx,sy,sz)

! cells reserved for the two port device
       
! ****** PORT OFFSETS TO BE GENERALISED -USE TERMINAL COORDINATES PERHAPS *******
       ixp(:)=ix
       iyp(:)=iy
       izp(:)=iz
              
       ixp(1)=ixp(1)-2*dx
       iyp(1)=iyp(1)-2*dy
       ixp(2)=ixp(2)-1*dx
       iyp(2)=iyp(2)-1*dy
       ixp(4)=ixp(4)+1*dx
       iyp(4)=iyp(4)+1*dy
       ixp(5)=ixp(5)+2*dx
       iyp(5)=iyp(5)+2*dy
       
       write(*,*)'Device cells:'
       do ii=1,5
         write(*,'(3I6)')ixp(ii),iyp(ii),izp(ii)
       end do
     
! cells of two port device cell centres
      write(*,*)'Circuit element:',i
      write(*,*)'Device positions and cells:'

      do ii=1,5
        CALL get_TLM_cell_centre_coordinate(ixp(ii),iyp(ii),izp(ii),sxp(ii),syp(ii),szp(ii))
        write(*,'(3ES16.6,3I6)')sxp(ii),syp(ii),szp(ii),ixp(ii),iyp(ii),izp(ii)
      end do

! Add the active elements, two TLM face patches to the surface list

      do ii=1,2
! Element 1 
        n_surfaces=n_surfaces+1
        surface_type(n_surfaces)=surface_type_zplane
        
        if (ii.EQ.1) then
          cell=2
        else
          cell=4
        end if
        
        write(*,*)'Setting active patch surface at cell:'
        write(*,*)ixp(cell),iyp(cell),izp(cell)
        
        surface_parameters(n_surfaces,1)=sxp(cell)-dl/2d0
        surface_parameters(n_surfaces,2)=syp(cell)-dl/2d0
        surface_parameters(n_surfaces,3)=szp(cell)-dl/2d0  ! zmin face
        surface_parameters(n_surfaces,4)=sxp(cell)+dl/2d0
        surface_parameters(n_surfaces,5)=syp(cell)+dl/2d0
        surface_parameters(n_surfaces,6)=szp(cell)-dl/2d0  ! zmin face

! Add the active element material of type 'SPICE' to the surface material list

        n_surface_materials=n_surface_materials+1  
    
        surface_material_type(n_surface_materials)=surface_material_type_SPICE  
    
        surface_material_to_surface_list(n_surface_materials)=n_surfaces
        SPICE_node_list(n_surface_materials)=ngspice_node_list(i,ii)      
        if (dy.NE.0) then
          SPICE_port_direction_list(n_surface_materials)='-y'        ! ****** PORT DIRECTION TO BE GENERALISED ******
        else
          SPICE_port_direction_list(n_surface_materials)='-x'        ! ****** PORT DIRECTION TO BE GENERALISED ******
        end if  
      
      end do ! next active patch
      
! Set the cells from each terminal to the active element cell to PEC
      do ii=1,ngspice_n_terminals(i)
      
        if (ii.EQ.1) then
          cell=1
        else if (ii.EQ.2) then
          cell=3
        else if (ii.eq.3) then
          cell=5
        end if
    
        write(*,*)'Add terminal connection cells:'
        write(*,'(3es16.6)')ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),ngspice_terminal_list(i,ii,3)
        write(*,'(3es16.6)')sxp(cell),syp(cell),szp(cell)
    
        CALL set_terminal_connection_cells(ngspice_terminal_list(i,ii,1),ngspice_terminal_list(i,ii,2),  &
                                           ngspice_terminal_list(i,ii,3),sxp(cell),syp(cell),szp(cell))
          
      end do
  
! Set the active element cells to free space - the active elements are set elsewhere
      material_mesh(centre,ixp(2),iyp(2),izp(2))=0
      material_mesh(centre,ixp(4),iyp(4),izp(4))=0
  
    else
    
      write(*,'(A54,I2,A10)')'ERROR: We do not yet hava a process for a device with ',ngspice_n_terminals(i),' terminals'
      STOP 1
  
    end if ! Number of terminals
    
! Read the package type
    
    read(10,'(A)')line
    
    if ( index(line,'none').NE.0 ) then
      
      package_type(i)=package_type_none
      
    else if ( index(line,'rec').NE.0 ) then
      
      package_type(i)=package_type_rectangular_block

! read the parameters
      
      read(10,*)(package_parameters(i,ii),ii=1,6)
      read(10,'(A)')material_name
! Add the dielectric block to the volume geometry list

      n_volumes=n_volumes+1
    
      if(n_volumes.GT.max_volumes) then
    
        write(*,*)'ERROR in GGI_TLM_create_PCB_simulation_model: maximum number of volumes exceeded'
        write(*,*)'Maximum number of volumes is set to ',max_volumes
        write(*,*)'in /GGI_TLM/SRC/TLM_MODULES/PCB_simulation_setup_modules.F90'
        STOP 1
    
      end if

      volume_type(n_volumes)=volume_type_rectangular_block2
    
      volume_parameters(n_volumes,1)=package_parameters(i,1)+sx
      volume_parameters(n_volumes,2)=package_parameters(i,2)+sy
      volume_parameters(n_volumes,3)=package_parameters(i,3)+sz
      volume_parameters(n_volumes,4)=package_parameters(i,4)+sx
      volume_parameters(n_volumes,5)=package_parameters(i,5)+sy
      volume_parameters(n_volumes,6)=package_parameters(i,6)+sz
    
! Volume material stuff

      n_volume_materials=n_volume_materials+1
      volume_material_type(n_volume_materials)=volume_material_type_DISPERSIVE
      volume_material_name(n_volume_materials)=trim(material_name)
      volume_material_to_volume_list(n_volume_materials)=n_volumes         
      
    else if ( index(line,'cyl').NE.0 ) then
      
      package_type(i)=package_type_cylinder
      
    else
    
      write(*,*)'ERROR: Unknown package type:',trim(line)
      STOP 1
    
    end if
  
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
