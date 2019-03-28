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
!     PROGRAM GGI_TLM_create_PCB_simulation_model
!
! DESCRIPTION
!     Convert template files for the simulation of a populted PCB in GGI_TLM
!     The PCB geometry is specified in gerber files
!     and components are added to this.
!     The electrical model of the components is solved using ngspice
!     and the electrical connection between comoponents is provied by GGI_TLM which is closely coupled to the ngspice solution.
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

PROGRAM GGI_TLM_create_PCB_simulation_model

USE TLM_general
USE File_information
USE PCB_simulation

IMPLICIT NONE

character(LEN=256) :: ipfilename
character(LEN=256) :: opfilename

! START

  CALL write_progress('STARTED: GGI_TLM_create_PCB_simulation_model')
  
!0. Reset variables

  n_gerber_files=0
  n_surfaces=0
  n_surface_materials=0
  n_PEC_surfaces=0
  n_volumes=0
  n_volume_materials=0
  tot_n_ngspice_nodes=0
  tot_n_ngspice_ports=0
  
!1a. Open the file containing the PCB simulation specifications

  write(*,*)'Enter the input filename for the PCB simulation specifications'
  read(*,'(A)')ipfilename
  
  OPEN(unit=10,file=trim(ipfilename),status='OLD',ERR=9000)

!1b. Open the template GGI_TLM input file

  read(10,'(A)')opfilename
  
  OPEN(unit=20,file=trim(opfilename))

!1c. Open the template Spice_circuit_TEMPLATE.cir
  
  OPEN(unit=30,file='Spice_circuit_TEMPLATE.cir')

!2  Process the PCB layout (gerber) files
!2a. Read gerber file and process to give stl file for GGI_TLM import
!2b. Read the z position for the layer
!2c. Add the layer to the Surface_list
!2d. Add the layer to the Sufrace_material_list as PEC to the Surface_list
!2d. Add the layer to the Sufrace_material_list as PEC

  CALL process_gerber()

!3. Add dielectric layers as required

!3a. Specift dielectric (xmin,ymin,zmin), (xmax,ymax,zmax) position in space and add to the Volume_list 
!3b. Choose a volume material file and add dielectric layers to the Volume_material_list 

  CALL add_dielectric_layers()

!4. Add vias as required
!4a. Read each via position (x,y,zmin,zmax) and add to the Surface_list
!4b. Add vias to the Surface_material_list as PEC
  
  ALLOCATE( material_mesh(1:7,1:mesh_nx+1,1:mesh_ny+1,1:mesh_nz+1) )  ! add an extra node on each dimension
  
  material_mesh(:,:,:,:)=0  ! reset the mesh

  CALL add_vias()

!5. Add components as required
!5a. Read component type 
!5b. Read number of ports  and the corresponding ngspice nodes
!5c. Read the component electrical connection positions in space (x1,y1,z1), (x2,y2,z2), ...
!5c. Generate the component geometry type and parameters
!5d. Write the component surfaces to the Surface_list 
!5e. Write the component volumes to the Volume_list
!5f. Write the component connections to the Surface_material_list and set as PEC
!5g. Write the component ngspice link surfaces to the Surface_material_list and set to type SPICE
!5h. Choose volume material files as required and add dielectric materials to the Volume_material_list 

  CALL add_components()

!6. Specify additional components 
!6a. Specify heatsink geometries as required and add to the Surface_list
!6b. Write the heatsink surfaces to the Surface_material_list and set as PEC

  CALL specify_additional_components()
  
!7 write the stl file for all PEC elements: vias, connections from PCB to Ngpsice link ports and additional components

  CALL mesh_to_stl()

!8 write the template GGI_TLM input file

  CALL write_GGI_TLM_input_file()

!9 write the template Spice_circuit_TEMPLATE.cir file

  CALL write_Spice_input_file()

!10 close files and finish 

  CLOSE(UNIT=10)
  
  CLOSE(UNIT=20)
  
  CLOSE(UNIT=30)
  
  DEALLOCATE( material_mesh )
  
  CALL write_progress('FINISHED: GGI_TLM_create_PCB_simulation_model')

  STOP
  
9000 write(*,*)'ERROR: Cannot open the file:',trim(ipfilename)
  STOP

END PROGRAM GGI_TLM_create_PCB_simulation_model
