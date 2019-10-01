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
! SUBROUTINE read_pml
!
! NAME
!     read_pml
!
! DESCRIPTION
!     read pml packet
!
! Example packet:
!
!PML
!4.0 4.0   4.0 4.0   4.0 4.0 ! PML thickness on each outer boundary surface: xmin,xmax,ymin,ymax,zmin,zmax
!1e-4                        ! PML reflection_coefficient
!
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!     1/11/2018 CJS:  Allow volumes to be created by filling surfaces
!
!
SUBROUTINE read_pml

USE TLM_general
USE file_information
USE geometry
USE PML_module
USE constants

IMPLICIT NONE

! local variables

! START  

  CALL write_line('CALLED: read_pml',0,output_to_screen_flag)
  
  PML_flag=.TRUE.
  
  read(input_file_unit,*,err=9000)pml_txmin,pml_txmax,pml_tymin,pml_tymax,pml_tzmin,pml_tzmax
      
  read(input_file_unit,*,err=9005)pml_r
  
  n_pml_volumes=0
  
  pml_volume_to_face(1:6)=0
  
  if (pml_txmin.GT.0.0) then
    n_pml_volumes=n_pml_volumes+1
    pml_volume_to_face(n_pml_volumes)=pml_face_xmin
    pml_xmin_flag=.TRUE.
  end if
  if (pml_txmax.GT.0.0) then
    n_pml_volumes=n_pml_volumes+1
    pml_volume_to_face(n_pml_volumes)=pml_face_xmax
    pml_xmax_flag=.TRUE.
  end if
  
  if (pml_tymin.GT.0.0) then
    n_pml_volumes=n_pml_volumes+1
    pml_volume_to_face(n_pml_volumes)=pml_face_ymin
    pml_ymin_flag=.TRUE.
  end if
  if (pml_tymax.GT.0.0) then
    n_pml_volumes=n_pml_volumes+1
    pml_volume_to_face(n_pml_volumes)=pml_face_ymax
    pml_ymax_flag=.TRUE.
  end if
  
  if (pml_tzmin.GT.0.0) then
    n_pml_volumes=n_pml_volumes+1
    pml_volume_to_face(n_pml_volumes)=pml_face_zmin
    pml_zmin_flag=.TRUE.
  end if
  if (pml_tzmax.GT.0.0) then
    n_pml_volumes=n_pml_volumes+1
    pml_volume_to_face(n_pml_volumes)=pml_face_zmax
    pml_zmax_flag=.TRUE.
  end if
  
  CALL write_line_integer('number of PML volumes',n_pml_volumes,0,output_to_screen_flag)

  CALL write_line('FINISHED: read_pml',0,output_to_screen_flag)
  
  RETURN
    
9000 CALL write_line('Error reading PML thicknesses from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
    
9005 CALL write_line('Error reading PML reflection parameter from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
END SUBROUTINE read_pml
