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
!SUBROUTINE get_cell_centre_coordinate(cell,point)
!SUBROUTINE get_cell_point_coordinate(face,point)
!SUBROUTINE get_cell_face_corner_coordinates(face,point1,point2,point3,point4)
!SUBROUTINE get_cell_corner_coordinates(cell,point1,point2,point3,point4,point5,point6,point7,point8)
!SUBROUTINE point_to_cell(point,cell)
!SUBROUTINE point_to_face(point,cell_point)
!SUBROUTINE point_to_cell_no_checks(point,cell)
!SUBROUTINE get_other_side_of_face(face1,face2)
!SUBROUTINE get_min_face(face1,face2)
!SUBROUTINE cell_point_inside_mesh(point,inside)
!FUNCTION   same_cell_point(point1,point2)
!SUBROUTINE get_cell_segment_face
!SUBROUTINE mesh_partition
!
! Name get_cell_centre_coordinate
!     
!
! Description
!     given the i,j,k cell number return the cell centre x,y,z coordinate
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_cell_centre_coordinate(cell,point)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(ijk)	:: cell
type(xyz)	:: point

! local_variables

! START

  point%x=mesh_xmin+cell%i*dl-dl/2d0
  point%y=mesh_ymin+cell%j*dl-dl/2d0
  point%z=mesh_zmin+cell%k*dl-dl/2d0

  RETURN
  
END SUBROUTINE get_cell_centre_coordinate
!
! Name get_cell_point_coordinate
!     
!
! Description
!     given the cell point number return the corresponding x,y,z coordinate
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_cell_point_coordinate(local_cell_point,xyz_point)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: local_cell_point
type(xyz)		:: xyz_point

! local_variables

real*8	:: offset_x
real*8	:: offset_y
real*8	:: offset_z

! START

  offset_x=-0.5d0
  offset_y=-0.5d0
  offset_z=-0.5d0
  
  if      (local_cell_point%point.eq.face_xmin) then
    offset_x=-1.0d0
  else if (local_cell_point%point.eq.face_xmax)  then
    offset_x= 0.0d0
  else if (local_cell_point%point.eq.face_ymin)  then
    offset_y=-1.0d0
  else if (local_cell_point%point.eq.face_ymax)  then
    offset_y= 0.0d0
  else if (local_cell_point%point.eq.face_zmin)  then
    offset_z=-1.0d0
  else if (local_cell_point%point.eq.face_zmax)  then
    offset_z= 0.0d0
  else if (local_cell_point%point.eq.centre)  then
    offset_x= -0.5d0
    offset_y= -0.5d0
    offset_z= -0.5d0
  else 
    write(*,*)'Error in get_local_cell_point_coordinates'
    write(*,*)'No face defined'
    write(*,*)'Face number (local_cell_point%point)=',local_cell_point%point
    STOP
  end if

  xyz_point%x=mesh_xmin+local_cell_point%cell%i*dl+dl*offset_x
  xyz_point%y=mesh_ymin+local_cell_point%cell%j*dl+dl*offset_y
  xyz_point%z=mesh_zmin+local_cell_point%cell%k*dl+dl*offset_z

  RETURN
  
END SUBROUTINE get_cell_point_coordinate
!
! Name get_cell_face_corner_coordinates
!     
!
! Description
!     given the face number return the 4 face corner x,y,z coordinates
!     note: the order of the points returned defines the normal direction
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_cell_face_corner_coordinates(face,point1,point2,point3,point4)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: face
type(xyz)		:: point1,point2,point3,point4

! local_variables

! START
 
  if (face%point.eq.face_xmin) then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point2%x=mesh_xmin+face%cell%i*dl-dl
    point2%y=mesh_ymin+face%cell%j*dl-dl
    point2%z=mesh_zmin+face%cell%k*dl
    point3%x=mesh_xmin+face%cell%i*dl-dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point4%x=mesh_xmin+face%cell%i*dl-dl
    point4%y=mesh_ymin+face%cell%j*dl
    point4%z=mesh_zmin+face%cell%k*dl-dl

  else if (face%point.eq.face_xmax)  then
  
    point1%x=mesh_xmin+face%cell%i*dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point4%x=mesh_xmin+face%cell%i*dl
    point4%y=mesh_ymin+face%cell%j*dl-dl
    point4%z=mesh_zmin+face%cell%k*dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point2%x=mesh_xmin+face%cell%i*dl
    point2%y=mesh_ymin+face%cell%j*dl
    point2%z=mesh_zmin+face%cell%k*dl-dl
    
  else if (face%point.eq.face_ymin)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point2%x=mesh_xmin+face%cell%i*dl
    point2%y=mesh_ymin+face%cell%j*dl-dl
    point2%z=mesh_zmin+face%cell%k*dl-dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl-dl
    point3%z=mesh_zmin+face%cell%k*dl
    point4%x=mesh_xmin+face%cell%i*dl-dl
    point4%y=mesh_ymin+face%cell%j*dl-dl
    point4%z=mesh_zmin+face%cell%k*dl
    
  else if (face%point.eq.face_ymax)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point4%x=mesh_xmin+face%cell%i*dl
    point4%y=mesh_ymin+face%cell%j*dl
    point4%z=mesh_zmin+face%cell%k*dl-dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point2%x=mesh_xmin+face%cell%i*dl-dl
    point2%y=mesh_ymin+face%cell%j*dl
    point2%z=mesh_zmin+face%cell%k*dl
    
  else if (face%point.eq.face_zmin)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point2%x=mesh_xmin+face%cell%i*dl-dl
    point2%y=mesh_ymin+face%cell%j*dl
    point2%z=mesh_zmin+face%cell%k*dl-dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl-dl
    point4%x=mesh_xmin+face%cell%i*dl
    point4%y=mesh_ymin+face%cell%j*dl-dl
    point4%z=mesh_zmin+face%cell%k*dl-dl
    
  else if (face%point.eq.face_zmax)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl
    point4%x=mesh_xmin+face%cell%i*dl-dl
    point4%y=mesh_ymin+face%cell%j*dl
    point4%z=mesh_zmin+face%cell%k*dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point2%x=mesh_xmin+face%cell%i*dl
    point2%y=mesh_ymin+face%cell%j*dl-dl
    point2%z=mesh_zmin+face%cell%k*dl
    
  else 
    write(*,*)'Error in get_cell_face_corner_coordinates'
    write(*,*)'No face defined'
    write(*,*)'face%point=',face%point
    STOP
  end if

  RETURN
  
END SUBROUTINE get_cell_face_corner_coordinates
!
! Name get_cell_corner_coordinates
!     
!
! Description
!     given the cell number return the 8 face corner x,y,z coordinates
!
! Comments:
!      
!
! History
!
!     started 30/08/12 CJS
!

SUBROUTINE get_cell_corner_coordinates(cell,point1,point2,point3,point4,point5,point6,point7,point8)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(ijk)	:: cell
type(xyz)	:: point1,point2,point3,point4,point5,point6,point7,point8

! local_variables

! START
 
    point1%x=mesh_xmin+cell%i*dl-dl
    point1%y=mesh_ymin+cell%j*dl-dl
    point1%z=mesh_zmin+cell%k*dl-dl
    point2%x=mesh_xmin+cell%i*dl-dl
    point2%y=mesh_ymin+cell%j*dl-dl
    point2%z=mesh_zmin+cell%k*dl
    point3%x=mesh_xmin+cell%i*dl-dl
    point3%y=mesh_ymin+cell%j*dl
    point3%z=mesh_zmin+cell%k*dl
    point4%x=mesh_xmin+cell%i*dl-dl
    point4%y=mesh_ymin+cell%j*dl
    point4%z=mesh_zmin+cell%k*dl-dl
  
    point5%x=mesh_xmin+cell%i*dl
    point5%y=mesh_ymin+cell%j*dl-dl
    point5%z=mesh_zmin+cell%k*dl-dl
    point6%x=mesh_xmin+cell%i*dl
    point6%y=mesh_ymin+cell%j*dl-dl
    point6%z=mesh_zmin+cell%k*dl
    point7%x=mesh_xmin+cell%i*dl
    point7%y=mesh_ymin+cell%j*dl
    point7%z=mesh_zmin+cell%k*dl
    point8%x=mesh_xmin+cell%i*dl
    point8%y=mesh_ymin+cell%j*dl
    point8%z=mesh_zmin+cell%k*dl-dl

  RETURN
  
END SUBROUTINE get_cell_corner_coordinates
!
! Name point_to_cell
!     
!
! Description
!      given the x,y,z coordinate return the i,j,k cell containing the point
!
! Comments:
!      
!
! History
!
!     started 10/08/12 CJS
!

SUBROUTINE point_to_cell(point,cell)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(xyz)	:: point
type(ijk)	:: cell

! local_variables

! START

  cell%i=int((point%x-mesh_xmin)/dl)+1
  cell%j=int((point%y-mesh_ymin)/dl)+1
  cell%k=int((point%z-mesh_zmin)/dl)+1
  
  if ((cell%i.lt.1).OR.(cell%i.gt.nx)) goto 9000
  if ((cell%j.lt.1).OR.(cell%j.gt.ny)) goto 9000
  if ((cell%k.lt.1).OR.(cell%k.gt.nz)) goto 9000
  
  RETURN

9000 write(*,*)'Error in point_to_cell'
     write(*,*)'The cell is outside the defined mesh'
     write(*,*)'x =',point%x,' y =',point%y,' z =',point%z
     write(*,*)'i =',cell%i ,' j =',cell%j ,' k =',cell%k
     write(*,*)'nx=',point%x,' ny=',point%y,' nz=',point%z
     STOP

  RETURN
  
END SUBROUTINE point_to_cell
!
! Name point_to_face
!     
!
! Description
!      given the x,y,z coordinate return the i,j,k cell enclosing the point and face closest to the point
!      in the cell_point structure
!
! Comments:
!      Used to determine the cell face for cable terminations
!
! History
!
!     started 19/11/12 CJS
!

SUBROUTINE point_to_face(point,cell_point_out)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(xyz)	:: point
type(cell_point):: cell_point_out

! local_variables

type(cell_point):: trial_cell_point
integer		:: face
type(xyz)	:: face_xyz
real*8		:: distance,min_distance

! function types

real*8		:: xyz_distance

! START

  cell_point_out%cell%i=int((point%x-mesh_xmin)/dl)+1
  cell_point_out%cell%j=int((point%y-mesh_ymin)/dl)+1
  cell_point_out%cell%k=int((point%z-mesh_zmin)/dl)+1
  
  if ((cell_point_out%cell%i.lt.1).OR.(cell_point_out%cell%i.gt.nx)) goto 9000
  if ((cell_point_out%cell%j.lt.1).OR.(cell_point_out%cell%j.gt.ny)) goto 9000
  if ((cell_point_out%cell%k.lt.1).OR.(cell_point_out%cell%k.gt.nz)) goto 9000

! loop over the 6 faces of the cell and find the closest face centre to the point

  trial_cell_point%cell=cell_point_out%cell
  min_distance=1e30
  
  do face=1,6  
  
    trial_cell_point%point=face
    
    CALL get_cell_point_coordinate(trial_cell_point,face_xyz)
  
    distance=xyz_distance(point,face_xyz)
    
    if (distance.lt.min_distance) then
      min_distance=distance
      cell_point_out=trial_cell_point
    end if
  
  end do
  
  RETURN

9000 write(*,*)'Error in point_to_face'
     write(*,*)'The cell is outside the defined mesh'
     write(*,*)'x =',point%x,' y =',point%y,' z =',point%z
     write(*,*)'i =',cell_point_out%cell%i ,' j =',cell_point_out%cell%j ,' k =',cell_point_out%cell%k
     write(*,*)'nx=',point%x,' ny=',point%y,' nz=',point%z
     STOP

  RETURN
  
END SUBROUTINE point_to_face
!
! Name point_to_cell_no_checks
!     
!
! Description
!      given the x,y,z coordinate return the i,j,k cell containing the point
!
! Comments:
!      
!
! History
!
!     started 10/08/12 CJS
!

SUBROUTINE point_to_cell_no_checks(point,cell)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(xyz)	:: point
type(ijk)	:: cell

! local_variables

! START

  cell%i=int((point%x-mesh_xmin)/dl)+1
  cell%j=int((point%y-mesh_ymin)/dl)+1
  cell%k=int((point%z-mesh_zmin)/dl)+1
  
  RETURN
  
END SUBROUTINE point_to_cell_no_checks
!
! Name get_other_side_of_face
!     
!
! Description
!     given the face number return the face centre x,y,z coordinate
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_other_side_of_face(face1,face2)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: face1
type(cell_point)	:: face2

! local_variables

! START

  face2%cell%i=face1%cell%i
  face2%cell%j=face1%cell%j
  face2%cell%k=face1%cell%k

  if      (face1%point.eq.face_xmin) then
    face2%cell%i=face1%cell%i-1
    face2%point=face_xmax
  else if (face1%point.eq.face_xmax)  then
    face2%cell%i=face1%cell%i+1
    face2%point=face_xmin
  else if (face1%point.eq.face_ymin)  then
    face2%cell%j=face1%cell%j-1
    face2%point=face_ymax 
  else if (face1%point.eq.face_ymax)  then
    face2%cell%j=face1%cell%j+1
    face2%point=face_ymin   
  else if (face1%point.eq.face_zmin)  then
    face2%cell%k=face1%cell%k-1
    face2%point=face_zmax 
  else if (face1%point.eq.face_zmax)  then
    face2%cell%k=face1%cell%k+1
    face2%point=face_zmin   
  else 
    write(*,*)'Error in get_other_side_of_faces'
    write(*,*)'No face defined'
    write(*,*)'Face number (face1%point)=',face1%point
    STOP
  end if

  RETURN
  
END SUBROUTINE get_other_side_of_face
!
! Name get_min_face
!     
!
! Description
!     given face1 return face2 on the min side of a cell (xmin,ymin or zmin as appropriate)
!
! Comments:
!      
!
! History
!
!     started 22/11/12 CJS
!

SUBROUTINE get_min_face(face1,face2)

USE geometry_types
USE cell_parameters

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: face1
type(cell_point)	:: face2

! local_variables

type(cell_point)	:: local_face

! START

! assume we have a min face to start with
  local_face%cell%i=face1%cell%i
  local_face%cell%j=face1%cell%j
  local_face%cell%k=face1%cell%k
  local_face%point=face1%point

  if      (face1%point.eq.face_xmin)  then
    local_face%cell%i=face1%cell%i
    local_face%point=face_xmin
  else if (face1%point.eq.face_ymin)  then
    local_face%cell%j=face1%cell%j
    local_face%point=face_ymin   
  else if (face1%point.eq.face_zmin)  then
    local_face%cell%k=face1%cell%k
    local_face%point=face_zmin   
  else if (face1%point.eq.face_xmax)  then
    local_face%cell%i=face1%cell%i+1
    local_face%point=face_xmin
  else if (face1%point.eq.face_ymax)  then
    local_face%cell%j=face1%cell%j+1
    local_face%point=face_ymin   
  else if (face1%point.eq.face_zmax)  then
    local_face%cell%k=face1%cell%k+1
    local_face%point=face_zmin   
  else 
    write(*,*)'Error in get_min_face'
    write(*,*)'cell=',local_face%cell%i,local_face%cell%j,local_face%cell%k
    write(*,*)'face number=',local_face%point
    if ( (local_face%point.GE.1).AND.(local_face%point.LE.7)) then
      write(*,*)'face=',face_string(local_face%point)
    end if
    STOP
  end if
  
  face2=local_face

  RETURN
  
END SUBROUTINE get_min_face
!
! Name cell_point_inside_mesh
!     
!
! Description
!     given a cell point, check that it is inside the mesh
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE cell_point_inside_mesh(local_cell_point,inside)

USE geometry_types
USE mesh
USE cell_parameters

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: local_cell_point
logical			:: inside

! local_variables

! START

  inside=.FALSE.

! first check, is the cell outside the defined mesh  
  if( (local_cell_point%cell%i.lt.1).OR.(local_cell_point%cell%i.gt.nx) ) RETURN
  if( (local_cell_point%cell%j.lt.1).OR.(local_cell_point%cell%j.gt.ny) ) RETURN
  if( (local_cell_point%cell%k.lt.1).OR.(local_cell_point%cell%k.gt.nz) ) RETURN
  
! second check, is the face on the outer boundary
  if( (local_cell_point%cell%i.eq.1).AND.(local_cell_point%point.eq.face_xmin) ) RETURN
  if( (local_cell_point%cell%j.eq.1).AND.(local_cell_point%point.eq.face_ymin) ) RETURN
  if( (local_cell_point%cell%k.eq.1).AND.(local_cell_point%point.eq.face_zmin) ) RETURN
  if( (local_cell_point%cell%i.eq.nx).AND.(local_cell_point%point.eq.face_xmax) ) RETURN
  if( (local_cell_point%cell%j.eq.ny).AND.(local_cell_point%point.eq.face_ymax) ) RETURN
  if( (local_cell_point%cell%k.eq.nz).AND.(local_cell_point%point.eq.face_zmax) ) RETURN
  
  inside=.TRUE.
  
  RETURN
  
END SUBROUTINE cell_point_inside_mesh
!
! Name same_cell_point
!     
!
! Description
!     check the equivalence of two cell_points.
!     Note the points are returned as equivalent if they are opposite sides of a face 
!
! Comments:
!      
!
! History
!
!     started 20/09/12 CJS
!

FUNCTION same_cell_point(cell_point1,cell_point2) RESULT(res)

USE geometry_types
USE cell_parameters

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: cell_point1
type(cell_point)	:: cell_point2
logical			:: res

! local_variables

! START

  res=.FALSE.

! first check, are the cells the same and the points the same? 
  if( (cell_point1%cell%i.EQ.cell_point2%cell%i).AND.	&
      (cell_point1%cell%j.EQ.cell_point2%cell%j).AND.	&
      (cell_point1%cell%k.EQ.cell_point2%cell%k).AND.	&
      (cell_point1%point.EQ.cell_point2%point) ) then      
    res=.TRUE.
    RETURN   
  end if
  
! check whether point1 and point2 are on opposite sides of a x directed face
  if( (cell_point1%cell%i+1.EQ.cell_point2%cell%i).AND.	&
      (cell_point1%cell%j.EQ.cell_point2%cell%j).AND.	&
      (cell_point1%cell%k.EQ.cell_point2%cell%k).AND.	&
      (cell_point1%point.EQ.face_xmax).AND.	   &
      (cell_point2%point.EQ.face_xmin) ) then	   
    res=.TRUE.
    RETURN   
  end if
  if( (cell_point1%cell%i-1.EQ.cell_point2%cell%i).AND.	&
      (cell_point1%cell%j.EQ.cell_point2%cell%j).AND.	&
      (cell_point1%cell%k.EQ.cell_point2%cell%k).AND.	&
      (cell_point1%point.EQ.face_xmin).AND.	   &
      (cell_point2%point.EQ.face_xmax) ) then	   
    res=.TRUE.
    RETURN   
  end if
  
! check whether point1 and point2 are on opposite sides of a y directed face
  if( (cell_point1%cell%i.EQ.cell_point2%cell%i).AND.	&
      (cell_point1%cell%j+1.EQ.cell_point2%cell%j).AND.	&
      (cell_point1%cell%k.EQ.cell_point2%cell%k).AND.	&
      (cell_point1%point.EQ.face_ymax).AND.	   &
      (cell_point2%point.EQ.face_ymin) ) then	   
    res=.TRUE.
    RETURN   
  end if
  if( (cell_point1%cell%i.EQ.cell_point2%cell%i).AND.	&
      (cell_point1%cell%j-1.EQ.cell_point2%cell%j).AND.	&
      (cell_point1%cell%k.EQ.cell_point2%cell%k).AND.	&
      (cell_point1%point.EQ.face_ymin).AND.	   &
      (cell_point2%point.EQ.face_ymax) ) then	   
    res=.TRUE.
    RETURN   
  end if
  
! check whether point1 and point2 are on opposite sides of a z directed face
  if( (cell_point1%cell%i.EQ.cell_point2%cell%i).AND.	&
      (cell_point1%cell%j.EQ.cell_point2%cell%j).AND.	&
      (cell_point1%cell%k+1.EQ.cell_point2%cell%k).AND.	&
      (cell_point1%point.EQ.face_zmax).AND.	   &
      (cell_point2%point.EQ.face_zmin) ) then	   
    res=.TRUE.
    RETURN   
  end if
  if( (cell_point1%cell%i.EQ.cell_point2%cell%i).AND.	&
      (cell_point1%cell%j.EQ.cell_point2%cell%j).AND.	&
      (cell_point1%cell%k-1.EQ.cell_point2%cell%k).AND.	&
      (cell_point1%point.EQ.face_zmin).AND.	   &
      (cell_point2%point.EQ.face_zmax) ) then	   
    res=.TRUE.
    RETURN   
  end if

  
END FUNCTION same_cell_point
!
! Name get_cell_segment_face
!     
!
! Description
!     given a cell segment, return the face number that the segment intersects
!
! Comments:
!      
!
! History
!
!     started 12/11/12 CJS
!

SUBROUTINE get_cell_segment_face(cell_segment_in,face)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(cell_segment)	:: cell_segment_in
integer			:: face

! local_variables

integer			:: face_point

! START

  if (cell_segment_in%segment_point(1)%point.eq.centre) then
! face point must be point 2. 
    face_point=2
  else
!  face point is point 1
    face_point=1
  end if
  
  face=cell_segment_in%segment_point(face_point)%point

! check we have a valid face  
  if ( (face.lt.1).OR.(face.gt.6) ) then
    write(*,*)'Error in get_cell_segment_face'
    write(*,*)'No face point found, face=',face
    write(*,*)'Cell segment point 1:'
    write(*,*)' i=',cell_segment_in%segment_point(1)%cell%i,	&
              ' j=',cell_segment_in%segment_point(1)%cell%j,	&
              ' k=',cell_segment_in%segment_point(1)%cell%k,	&
              ' face=',cell_segment_in%segment_point(1)%point
    write(*,*)'Cell segment point 2:'
    write(*,*)' i=',cell_segment_in%segment_point(2)%cell%i,	&
              ' j=',cell_segment_in%segment_point(2)%cell%j,	&
              ' k=',cell_segment_in%segment_point(2)%cell%k,	&
              ' face=',cell_segment_in%segment_point(2)%point
    STOP
  end if
  
  RETURN
  
END SUBROUTINE get_cell_segment_face
!
! Name mesh_partition
!     
!
! Description
!   partition mesh for parallel solver - divides up the mesh in the z direction     
!
! Comments:
!      
!
! History
!
!     started 9/07/09 CJS
!

SUBROUTINE mesh_partition( )

USE Mesh
USE TLM_general
USE File_information

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer :: cell,p_rank,last_p_rank
  integer :: face

! function_types

! START

  CALL write_line('CALLED: mesh_partition',0,output_to_screen_flag)

  ALLOCATE(cell_rank(1:nz))

  ALLOCATE(cell_face_rank(1:nz,1:6))
  
! work out the z extent of this mesh partition
 
  if (np.gt.nz) goto 9000

! set the cell ranks of every cell in the z direction

  do cell=1,nz
  
    cell_rank(cell)=INT(np*real(cell-1)/real(nz))
    
  end do
  
  nz1=-1
  nz2=-1
  
  do cell=1,nz
  
    if ( (cell_rank(cell).eq.rank).AND.nz1.eq.-1) then
      nz1=cell
    end if
    if (  cell_rank(cell).eq.rank)then
      nz2=cell
    end if
  
  end do
  
! allocate the mesh 1 cell beyond the range owned by this processor 
! unless we are on the outer boundary

  nzmin=nz1-1  
  nzmin=max(1,nzmin)
  
  nzmax=nz2+1  
  nzmax=min(nz,nzmax)

! OLD
!! this process should own the zmin face of the nzmax cell and the zmax face of the nzmin cell
!  do cell=1,nz  
!    if (cell.eq.nzmax) then
!      cell_face_rank(cell,face_zmin)=rank
!    end if   
!    if (cell.eq.nzmin) then
!      cell_face_rank(cell,face_zmax)=rank
!    end if
!  end do

! this process should own the zmin face of the nzmax cell and the zmax face of the nzmin cell
! set cell_face_rank to cell_rank initially    
  do cell=1,nz  
    do face=1,6
      cell_face_rank(cell,face)=cell_rank(cell)
    end do
  end do
   
  do cell=1,nz-1
  
! look for a change in the rank of the process across cells in z

    if (cell_rank(cell).NE.cell_rank(cell+1)) then
! we are at the interface between processes, here the zmax face of the cell belongs to cell+1

      cell_face_rank(cell,face_zmax)=cell_rank(cell+1)
      
    end if
    
  end do ! next cell
  
  write(*,*)'Rank=',rank,' np=',np,' nzmin=',nzmin,' nz1=',nz1,'nz2=',nz2,' nzmax=',nzmax
  
  CALL write_line('FINISHED: mesh_partition',0,output_to_screen_flag)
 
  RETURN

9000 write(*,*)'Error in mesh_partition'
     write(*,*)'Number of processors greater than z dimension'
     write(*,*)'Np=',np,' nz=',nz

     STOP
  
END SUBROUTINE mesh_partition
