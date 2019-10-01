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
!SUBROUTINE build_PML_geometry
!SUBROUTINE build_PML_volume
!SUBROUTINE plot_PML_tet_volumes
!
! NAME
!     SUBROUTINE build_PML_geometry
!
! DESCRIPTION
!     build_PML_geometry:
!
!     Create volume geometric entities for PML volumes as required
!     The geometric entities consist of trinagulated volumes
!     Triangulated volumes are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started  1/10/2019 CJS  Add PML
!
!
!
SUBROUTINE build_PML_geometry

USE TLM_general
USE geometry_types
USE geometry
USE PML_module
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number

! START

  CALL write_line('CALLED: build_PML_geometry',0,output_to_screen_flag)

! Allocate the volumes for the PML geometry

  ALLOCATE( pml_volumes(1:n_pml_volumes) )

  volume_number=0
  
  if (pml_xmin_flag) then  
    volume_number=volume_number+1
    CALL build_PML_volume(volume_number,mesh_xmin,mesh_xmin+pml_txmin,mesh_ymin,mesh_ymax,mesh_zmin,mesh_zmax)    
  end if

  if (pml_xmax_flag) then  
    volume_number=volume_number+1
    CALL build_PML_volume(volume_number,mesh_xmax-pml_txmax,mesh_xmax,mesh_ymin,mesh_ymax,mesh_zmin,mesh_zmax)    
  end if

  if (pml_ymin_flag) then  
    volume_number=volume_number+1
    CALL build_PML_volume(volume_number,mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymin+pml_tymin,mesh_zmin,mesh_zmax)    
  end if

  if (pml_ymax_flag) then  
    volume_number=volume_number+1
    CALL build_PML_volume(volume_number,mesh_xmin,mesh_xmax,mesh_ymax-pml_tymax,mesh_ymax,mesh_zmin,mesh_zmax)    
  end if

  if (pml_zmin_flag) then  
    volume_number=volume_number+1
    CALL build_PML_volume(volume_number,mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymax,mesh_zmin,mesh_zmin+pml_tzmin)    
  end if

  if (pml_zmax_flag) then  
    volume_number=volume_number+1
    CALL build_PML_volume(volume_number,mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymax,mesh_zmax-pml_tzmax,mesh_zmax)    
  end if
  
  if (write_geometry_vtk_files) then
    CALL plot_PML_tet_volumes()
  end if
  
  CALL write_line('FINISHED: build_PML_geometry',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error in build_PML_geometry:',0,.TRUE.)
     CALL write_line('volume type number not defined',0,.TRUE.)
     CALL write_line_integer('volume type number',pml_volumes(volume_number)%volume_type,0,.TRUE.)
     STOP

  
END SUBROUTINE build_PML_geometry
!
! NAME
!     SUBROUTINE build_PML_volume
!
! DESCRIPTION
!     build_PML_volume
!
!     create a triangulated rectangular block volume for a PML volume
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/10/2019 CJS  Add PML
!
!
SUBROUTINE build_PML_volume(volume_number,x1in,x2in,y1in,y2in,z1in,z2in)


USE TLM_general
USE geometry_types
USE geometry
USE PML_module
USE file_information
USE constants

IMPLICIT NONE

integer	:: volume_number
real*8	:: x1in,y1in,z1in,x2in,y2in,z2in

! local variables

integer	:: number_of_tets
integer	:: tet_count

real*8	:: x1,y1,z1,x2,y2,z2
real*8	:: swap

type(xyz)	:: point0
type(xyz)	:: point1,point2,point3,point4
type(xyz)	:: point5,point6,point7,point8

! START

! CREATE A TRIANGULATED RECTANGULAR BLOCK DEFINED BY OPPOSITE CORNER COORDINATES
    
      number_of_tets=12
      pml_volumes(volume_number)%number_of_tets=number_of_tets    
      allocate( pml_volumes(volume_number)%tet_list(1:number_of_tets) )
      
      pml_volumes(volume_number)%volume_type=volume_type_rectangular_block2
      
      x1=x1in
      y1=y1in
      z1=z1in
      x2=x2in
      y2=y2in
      z2=z2in
      
! rectangular_block vertices
      pml_volumes(volume_number)%volume_parameters(1)=x1
      pml_volumes(volume_number)%volume_parameters(2)=y1
      pml_volumes(volume_number)%volume_parameters(3)=z1
      pml_volumes(volume_number)%volume_parameters(4)=x2
      pml_volumes(volume_number)%volume_parameters(5)=y2
      pml_volumes(volume_number)%volume_parameters(6)=z2

! transformation      
      pml_volumes(volume_number)%trans%parameters(1)=0d0
      pml_volumes(volume_number)%trans%parameters(2)=0d0
      pml_volumes(volume_number)%trans%parameters(3)=0d0
      pml_volumes(volume_number)%trans%parameters(4)=0d0
      pml_volumes(volume_number)%trans%parameters(5)=0d0
      pml_volumes(volume_number)%trans%parameters(6)=0d0
      
! ensure that x2>x1, y2>y1 and z2>z1. This ensures that the volume normals are outward

      if (x1.GT.x2) then
        swap=x1
	x1=x2
	x2=swap
      end if
      if (y1.GT.y2) then
        swap=y1
	y1=y2
	y2=swap
      end if
      if (z1.GT.z2) then
        swap=z1
	z1=z2
	z2=swap
      end if

! create central point, point0

      point0%x=(x1+x2)/2d0
      point0%y=(y1+y2)/2d0
      point0%z=(z1+z2)/2d0

      point1%x=x1
      point1%y=y1
      point1%z=z1
      
      point2%x=x2
      point2%y=y1
      point2%z=z1
      
      point3%x=x2
      point3%y=y2
      point3%z=z1
      
      point4%x=x1
      point4%y=y2
      point4%z=z1
      
      point5%x=x1
      point5%y=y1
      point5%z=z2
      
      point6%x=x2
      point6%y=y1
      point6%z=z2
      
      point7%x=x2
      point7%y=y2
      point7%z=z2
      
      point8%x=x1
      point8%y=y2
      point8%z=z2
    	  
! create volume tets    
! set tets to make the normal point outwards	  
      tet_count=0
	  
      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point3
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point2
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point4
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point3
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point2
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point6
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point6
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point5
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point3
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point7
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point2
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point7
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point6
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point5
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point6
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point7
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point5
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point7
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point8
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point4
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point7
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point3
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point4
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point8
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point7
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0
	  
      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point8
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point4
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

      tet_count=tet_count+1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(1)=point1
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(2)=point5
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(3)=point8
      pml_volumes(volume_number)%tet_list(tet_count)%vertex(4)=point0

  RETURN

END SUBROUTINE build_PML_volume
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
! NAME
!     SUBROUTINE plot_PML_tet_volumes
!
! DESCRIPTION
!     plot_PML_tet_volumes:
!
!     All the tet volumes are written to .vtk files for visualisation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE plot_PML_tet_volumes()

USE TLM_general
USE geometry_types
USE geometry
USE PML_module
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number

integer	:: number_of_tets

! START

  CALL write_line('CALLED: plot_PML_tet_volumes',0,output_to_screen_flag)

  do volume_number=1,n_pml_volumes

    number_of_tets=pml_volumes(volume_number)%number_of_tets

    if (number_of_tets.GT.0) then

! open and write triangulated volume to vtk format file
      CALL open_vtk_file(tet_volume_file_unit,pml_tet_volume_file_extension,volume_number) 
      
      CALL write_volume_list_vtk(tet_volume_file_unit,	&
                                  number_of_tets,pml_volumes(volume_number)%tet_list)
      
      CALL close_vtk_file(tet_volume_file_unit) 
    
    end if ! number_of_tets.gt.0

  end do ! next volume number

  CALL write_line('FINISHED: plot_PML_tet_volumes',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE plot_PML_tet_volumes
