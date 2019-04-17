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
!     SUBROUTINE process_gerber
!
! DESCRIPTION
!2  Process the PCB layout (gerber) files
!2a. Read gerber file and process to give stl file for GGI_TLM import
!2b. Read the z position for the layer
!2c. Add the layer
!2a. Read gerber file and process to give stl file for GGI_TLM import
!2b. Read the z position for the layer
!2c. Add the layer to the Surface_list
!2d. Add the layer to the Sufrace_material_list as PEC to the Surface_list
!2d. Add the layer to the Sufrace_material_list as PEC
!	
!     
! COMMENTS
!     
! Example of the relevant section of the input file:
!
! 1      # Number of gerber files to include
! two_track.gbr
! 0.0    # z position for the layer specified in this gerber file
!
!
!
! HISTORY
!
!     started 14/3/19 CJS
!     
!

SUBROUTINE process_gerber

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i

character(LEN=256) :: command
character(LEN=16) :: number

integer :: iout
logical :: file_exists

! START

!2  Process the PCB layout (gerber) files

read(10,*)dl

read(10,*) xmin,ymin,zmin,xmax,ymax,zmax     ! TLM problem space dimensions

! calculate the number of cells in each direction
mesh_nx=NINT((xmax-xmin)/dl)
mesh_ny=NINT((ymax-ymin)/dl)
mesh_nz=NINT((zmax-zmin)/dl)

! adjust max dimensions to ensure there is an integer number of cells in each direction
xmax=xmin+mesh_nx*dl
ymax=ymin+mesh_ny*dl
zmax=zmin+mesh_nz*dl

read(10,*)n_gerber_files
write(*,*)'Number of gerber files=',n_gerber_files

if (n_gerber_files.GT.max_gerber_files) then
  write(*,*)'ERROR in GGI_TLM_create_PCB_simulation_model: maximum number of gerber files exceeded'
  write(*,*)'Maximum number of gerber files is set to ',max_gerber_files
  write(*,*)'in /GGI_TLM/SRC/TLM_MODULES/PCB_simulation_setup_modules.F90'
  STOP
end if

do i=1,n_gerber_files

!2a. Read gerber file and process to give stl file for GGI_TLM import
  read(10,'(A)')gerber_filename(i)
  
! Check that the file exists

  INQUIRE(FILE=trim(gerber_filename(i)), EXIST=file_exists)
  
  if (.NOT.file_exists) GOTO 9000
  
  write(number,'(ES16.6)')dl/4d0   ! gerber to stl resolution - should be less than dl
  command=trim(executable_path)//'GGI_TLM_gerber_to_stl '//trim(gerber_filename(i))//' '//trim(number)
  
  write(*,*)'Running gerber to stl command:',trim(command)
  CALL execute_command_line (trim(command), exitstat=iout)
  print *, "Exit status of GGI_TLM_gerber_to_stl was ", iout 
  
  n_surfaces=n_surfaces+1

!2b. Read the z position for the layer
  read(10,*)gerber_z_offset(i)
  surface_z_offset(n_surfaces)=gerber_z_offset(i)
  
!2c. Add the layer to the Surface_list

  surface_type(n_surfaces)=surface_type_stl
  surface_filename(n_surfaces)=trim(gerber_filename(i))//'.stl'

!2d. Add the layer to the Sufrace_material_list as PEC to the Surface_list

  if (n_PEC_surfaces.EQ.0) then
    n_surface_materials=n_surface_materials+1  
    surface_material_type(n_surface_materials)=surface_material_type_PEC
  end if
  n_PEC_surfaces=n_PEC_surfaces+1
  PEC_surface_list(n_PEC_surfaces)=n_surfaces
  PEC_surface_orientation_list(n_PEC_surfaces)=1
end do

RETURN  
  
9000 write(*,*)'ERROR: Cannot open the file:',trim(gerber_filename(i))
  STOP

END SUBROUTINE process_gerber
