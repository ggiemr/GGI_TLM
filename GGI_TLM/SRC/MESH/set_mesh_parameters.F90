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
! SUBROUTINE set_mesh_parameters
!
! NAME
!     set_mesh_parameters
!
! DESCRIPTION
!     
!     Set the general mesh parameters:
!     1. mesh dimension in cells, nx,ny,nz
!     2. mesh limits mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymax,mesh_zmin,mesh_zmax
!     3. cell dimension, dl
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE set_mesh_parameters()

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

  real*8 dlx,dly,dlz

! START  
  
  CALL write_line('CALLED: set_mesh_parameters',0,output_to_screen_flag)
  
  if (dl.eq.0d0) then
  
! we have not set the mesh edge length so must have mesh dimension in cells defined 

    if ( (nx.eq.0).OR.(ny.eq.0).OR.(nz.eq.0) ) GOTO 9000

! calculate dl
    dlx=(mesh_xmax-mesh_xmin)/dble(nx)
    dly=(mesh_ymax-mesh_ymin)/dble(ny)
    dlz=(mesh_zmax-mesh_zmin)/dble(nz)

! calculate dl as the average of dlx, dly, dlz    
    dl=(dlx+dly+dlz)/(3d0)

! write warning if the defined mesh is not cubic    
    if ( (dlx.ne.dl).OR.(dly.ne.dl).OR.(dlz.ne.dl) ) then
    
      write(warning_file_unit,*)'Set_mesh_parameters:'
      write(warning_file_unit,*)'(dlx.ne.dl).OR.(dly.ne.dl).OR.(dlz.ne.dl)'
      write(warning_file_unit,*)'dl=',dl
      write(warning_file_unit,*)'dlx=',dlx
      write(warning_file_unit,*)'dly=',dly
      write(warning_file_unit,*)'dlz=',dlz
      
    end if

! reset nx,ny,nz for a uniform grid
    nx=nint((mesh_xmax-mesh_xmin)/dl)
    ny=nint((mesh_ymax-mesh_ymin)/dl)
    nz=nint((mesh_zmax-mesh_zmin)/dl)
  
  else
! we have set dl so must work out nx,ny,nz

    if ( (mesh_xmin.eq.0d0).AND.(mesh_xmax.eq.0d0).AND.	&
         (mesh_ymin.eq.0d0).AND.(mesh_ymax.eq.0d0).AND.	&
         (mesh_zmin.eq.0d0).AND.(mesh_zmax.eq.0d0) ) GOTO 9010

    nx=nint((mesh_xmax-mesh_xmin)/dl)
    ny=nint((mesh_ymax-mesh_ymin)/dl)
    nz=nint((mesh_zmax-mesh_zmin)/dl)

! make sure that nx, ny and nz are greater than or equal to 1
    
    if (nx.eq.0) nx=1
    if (ny.eq.0) ny=1
    if (nz.eq.0) nz=1
    
  end if  
  
  CALL write_line('FINISHED: set_mesh_parameters',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error in set_mesh_parameters:',0,.TRUE.)
     CALL write_line('dl unset and nx,ny or nz=0',0,.TRUE.)
     STOP
  
9010 CALL write_line('Error in set_mesh_parameters:',0,.TRUE.)
     CALL write_line('dl set and mesh_xmin=mesh_xmax=mesh_ymin=mesh_ymax=mesh_zmin=mesh_zmax=0.0',0,.TRUE.)
     STOP

  
  RETURN
  
END SUBROUTINE set_mesh_parameters
