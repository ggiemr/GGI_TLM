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
! Example packet 1:
!
!PML
!4.0  4.0   4.0  4.0   4.0 4.0   ! PML thickness on each outer boundary surface: xmin,xmax,ymin,ymax,zmin,zmax
!1e-4 1e-4  1e-4 1e-4  1e-4 1e-4 ! PML reflection_coefficient on each outer boundary surface: xmin,xmax,ymin,ymax,zmin,zmax
!2                               ! PML order
!
! Example packet 2:
!
!PML
!4.0   ! PML thickness on all outer boundary surfaces
!1e-4  ! PML reflection_coefficient on all outer boundary surfaces
!2     ! PML order
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

  character(LEN=256) :: line

! START  

  CALL write_line('CALLED: read_pml',0,output_to_screen_flag)
  
  if (PML_flag) then
    write(*,*)'ERROR: already read PML information from the input file'
    STOP 1
  end if
  
  PML_flag=.TRUE.
  
  read(input_file_unit,'(A)',err=9000,end=9000)line
  
  ! attempt to read 6 values for PML layer thickness. If there are not six values, read one value
  ! and use it for all  mesh boundaries
  
  read(line,*,err=1000)pml_txmin,pml_txmax,pml_tymin,pml_tymax,pml_tzmin,pml_tzmax
  GOTO 1010 ! read six values OK
  
1000 CONTINUE
     read(line,*,err=9000,end=9000)pml_txmin
     pml_txmax=pml_txmin
     pml_tymin=pml_txmin
     pml_tymax=pml_txmin
     pml_tzmin=pml_txmin
     pml_tzmax=pml_txmin

1010 CONTINUE

! Checks on values
  if ( (pml_txmin.LT.0.0).OR.(pml_txmax.LT.0.0).OR.    &
       (pml_tymin.LT.0.0).OR.(pml_tymax.LT.0.0).OR.    &
       (pml_tzmin.LT.0.0).OR.(pml_tzmax.LT.0.0) ) then
       
    write(*,*)'ERROR in read_pml: PML layer thickness should be greater than or equal to zero'
    write(*,*)'pml_txmin=',pml_txmin
    write(*,*)'pml_txmax=',pml_txmax
    write(*,*)'pml_tymin=',pml_tymin
    write(*,*)'pml_tymax=',pml_tymax
    write(*,*)'pml_tzmin=',pml_tzmin
    write(*,*)'pml_tzmax=',pml_tzmax
    
    STOP 1
       
  end if
 
  read(input_file_unit,'(A)',err=9000,end=9000)line
  
  ! attempt to read 6 values for PML layer reflection coefficient. If there are not six values, read one value
  ! and use it for all  mesh boundaries
  
  read(line,*,err=2000)pml_r_xmin,pml_r_xmax,pml_r_ymin,pml_r_ymax,pml_r_zmin,pml_r_zmax
  GOTO 2010 ! read six values OK
  
2000 CONTINUE
     read(line,*,err=9005,end=9005)pml_r_xmin
     pml_r_xmax=pml_r_xmin
     pml_r_ymin=pml_r_xmin
     pml_r_ymax=pml_r_xmin
     pml_r_zmin=pml_r_xmin
     pml_r_zmax=pml_r_xmin

2010 CONTINUE

! Checks on values
  if ( (pml_r_xmin.LE.0.0).OR.(pml_r_xmax.LE.0.0).OR.    &
       (pml_r_ymin.LE.0.0).OR.(pml_r_ymax.LE.0.0).OR.    &
       (pml_r_zmin.LE.0.0).OR.(pml_r_zmax.LE.0.0).OR.    & 
       (pml_r_xmin.GE.1.0).OR.(pml_r_xmax.GE.1.0).OR.    &
       (pml_r_ymin.GE.1.0).OR.(pml_r_ymax.GE.1.0).OR.    &
       (pml_r_zmin.GE.1.0).OR.(pml_r_zmax.GE.1.0) ) then
       
    write(*,*)'ERROR in read_pml: PML reflectivity should be between 0 and 1'
    write(*,*)'pml_r_xmin=',pml_r_xmin
    write(*,*)'pml_r_xmax=',pml_r_xmax
    write(*,*)'pml_r_ymin=',pml_r_ymin
    write(*,*)'pml_r_ymax=',pml_r_ymax
    write(*,*)'pml_r_zmin=',pml_r_zmin
    write(*,*)'pml_r_zmax=',pml_r_zmax
    
    STOP 1
       
  end if
                
  read(input_file_unit,*,err=9010)pml_order
  
  if (pml_order.LE.0) then
       
    write(*,*)'ERROR in read_pml: PML order should be greater than 0'
    write(*,*)'pml_order=',pml_order
    
    STOP 1
       
  end if
  
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
     STOP 1
    
9005 CALL write_line('Error reading PML reflection parameter from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
      
9010 CALL write_line('Error reading PML order from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP 1
 
END SUBROUTINE read_pml
