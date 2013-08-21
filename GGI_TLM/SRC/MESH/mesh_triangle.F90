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
!SUBROUTINE mesh_triangle
!
! NAME
!     SUBROUTINE mesh_triangle
!
! DESCRIPTION
!     mesh_triangle:
!
!     mesh a single triangle
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/08/2012 CJS
!     revised 14/9/2012 CJS - simpler strategy which avoids duplicate faces
!
SUBROUTINE mesh_triangle(surface_mesh,triangle)
  
USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE file_information
USE constants

IMPLICIT NONE

  integer		:: surface_mesh(1:nx,1:ny,1:nz,1:6)
  type(xyz_triangle)	:: triangle

! local variables

  integer 		:: ix,iy,iz
   
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
 
  type(xyz)		:: min_point,max_point
  type(ijk)		:: min_cell,max_cell
  
  type(xyz)		:: point1,point2
  type(cell_point)	:: cell_point1,cell_point2
  
  type(cell_point)	:: mesh_face
  
  type(xyz)		:: intersection_point
  type(ijk)		:: intersection_cell
  logical		:: intersection
     
  type(xyz)		:: normal
  
  real*8		:: dmin,dmax
  
  logical		:: inside
 
! function types  
  real*8		:: xyz_distance
  
! START
  
  x1=triangle%vertex(1)%x
  y1=triangle%vertex(1)%y
  z1=triangle%vertex(1)%z
  
  x2=triangle%vertex(2)%x
  y2=triangle%vertex(2)%y
  z2=triangle%vertex(2)%z
  
  x3=triangle%vertex(3)%x
  y3=triangle%vertex(3)%y
  z3=triangle%vertex(3)%z

! get the extent of the triangle in x,y and z
  min_point%x=min(x1,x2,x3)
  min_point%y=min(y1,y2,y3)
  min_point%z=min(z1,z2,z3)
  
  max_point%x=max(x1,x2,x3)
  max_point%y=max(y1,y2,y3)
  max_point%z=max(z1,z2,z3)
    
! get the extent of the triangle in cells in x,y and z, note cells can be outside the mesh here
  CALL point_to_cell_no_checks(min_point,min_cell)
  CALL point_to_cell_no_checks(max_point,max_cell)
  
! LOOK FOR FACES NORMAL TO X
  do iy=min_cell%j,max_cell%j 
    do iz=min_cell%k,max_cell%k
! get two points on the min and max x extent of the mesh, through the centre of the x face of cell iy,iz
      
      cell_point1%cell%i=1
      cell_point1%cell%j=iy
      cell_point1%cell%k=iz
      cell_point1%point=face_xmin
      
      cell_point2%cell%i=nx
      cell_point2%cell%j=iy
      cell_point2%cell%k=iz
      cell_point2%point=face_xmax
      
      CALL get_cell_point_coordinate(cell_point1,point1)
                
      CALL get_cell_point_coordinate(cell_point2,point2)

! check whether the point defined by point1,point2 intersects the current triangle

      CALL line_triangle_intersection(point1,point2,triangle,	&
	                              intersection_point,intersection )
	       
      if (intersection) then
	
        CALL point_to_cell_no_checks(intersection_point,intersection_cell)
	
        mesh_face%cell%i=intersection_cell%i
        mesh_face%cell%j=intersection_cell%j
        mesh_face%cell%k=intersection_cell%k

! we must now find the closest normal to x face to the intersection point	
        cell_point1%cell%i=intersection_cell%i
        cell_point1%cell%j=intersection_cell%j
        cell_point1%cell%k=intersection_cell%k
        cell_point1%point=face_xmin
        CALL get_cell_point_coordinate(cell_point1,point1)
      
        cell_point2%cell%i=intersection_cell%i
        cell_point2%cell%j=intersection_cell%j
        cell_point2%cell%k=intersection_cell%k
        cell_point2%point=face_xmax
        CALL get_cell_point_coordinate(cell_point2,point2)
      
        dmin=xyz_distance(point1,intersection_point)
        dmax=xyz_distance(point2,intersection_point)
	
	if (dmin.lt.dmax) then
          mesh_face%point=face_xmin
	else
          mesh_face%point=face_xmax	
	end if
	
! calculate triangle normal and swap to the other side of the face if required to keep the mesh normal the same
        CALL triangle_normal(triangle,normal)
	if ( (normal%x.gt.0d0).AND.(mesh_face%point.EQ.face_xmin) ) then
          mesh_face%cell%i=mesh_face%cell%i-1
	  mesh_face%point=face_xmax
	else if ( (normal%x.lt.0d0).AND.(mesh_face%point.EQ.face_xmax) ) then
          mesh_face%cell%i=mesh_face%cell%i+1
	  mesh_face%point=face_xmin
	end if
	
	CALL cell_point_inside_mesh(mesh_face,inside)
	
	if (inside) then ! set face in the surface mesh array
	  surface_mesh(mesh_face%cell%i,mesh_face%cell%j,mesh_face%cell%k,mesh_face%point)=1
	end if
	
      end if
      
    end do ! next z cell
  end do ! next y cell 
  
! LOOK FOR FACES NORMAL TO Y
  do ix=min_cell%i,max_cell%i 
    do iz=min_cell%k,max_cell%k
! get two points on the min and max x extent of the mesh, through the centre of the x face of cell iy,iz
      
      cell_point1%cell%i=ix
      cell_point1%cell%j=1
      cell_point1%cell%k=iz
      cell_point1%point=face_ymin
      
      cell_point2%cell%i=ix
      cell_point2%cell%j=ny
      cell_point2%cell%k=iz
      cell_point2%point=face_ymax
      
      CALL get_cell_point_coordinate(cell_point1,point1)
      CALL get_cell_point_coordinate(cell_point2,point2)

! check whether the point defined by point1,point2 intersects the current triangle

      CALL line_triangle_intersection(point1,point2,triangle,	&
	                              intersection_point,intersection )
	       
      if (intersection) then
	
        CALL point_to_cell_no_checks(intersection_point,intersection_cell)
	
        mesh_face%cell%i=intersection_cell%i
        mesh_face%cell%j=intersection_cell%j
        mesh_face%cell%k=intersection_cell%k

! we must now find the closest normal to y face to the intersection point	
        cell_point1%cell%i=intersection_cell%i
        cell_point1%cell%j=intersection_cell%j
        cell_point1%cell%k=intersection_cell%k
        cell_point1%point=face_ymin
        CALL get_cell_point_coordinate(cell_point1,point1)
      
        cell_point2%cell%i=intersection_cell%i
        cell_point2%cell%j=intersection_cell%j
        cell_point2%cell%k=intersection_cell%k
        cell_point2%point=face_ymax
        CALL get_cell_point_coordinate(cell_point2,point2)
      
        dmin=xyz_distance(point1,intersection_point)
        dmax=xyz_distance(point2,intersection_point)
	
	if (dmin.lt.dmax) then
          mesh_face%point=face_ymin
	else
          mesh_face%point=face_ymax	
	end if
	
! calculate triangle normal and swap to the other side of the face if required to keep the mesh normal the same
        CALL triangle_normal(triangle,normal)
	if ( (normal%y.gt.0d0).AND.(mesh_face%point.EQ.face_ymin) ) then
          mesh_face%cell%j=mesh_face%cell%j-1
	  mesh_face%point=face_ymax
	else if ( (normal%y.lt.0d0).AND.(mesh_face%point.EQ.face_ymax) ) then
          mesh_face%cell%j=mesh_face%cell%j+1
	  mesh_face%point=face_ymin
	end if
	
	CALL cell_point_inside_mesh(mesh_face,inside)
	
	if (inside) then ! set face in the surface mesh array
	  surface_mesh(mesh_face%cell%i,mesh_face%cell%j,mesh_face%cell%k,mesh_face%point)=1
	end if
	
      end if
      
    end do ! next z cell
  end do ! next y cell 
  
! LOOK FOR FACES NORMAL TO Z
  do ix=min_cell%i,max_cell%i 
    do iy=min_cell%j,max_cell%j
! get two points on the min and max x extent of the mesh, through the centre of the x face of cell iy,iz
      
      cell_point1%cell%i=ix
      cell_point1%cell%j=iy
      cell_point1%cell%k=1
      cell_point1%point=face_zmin
      
      cell_point2%cell%i=ix
      cell_point2%cell%j=iy
      cell_point2%cell%k=nz
      cell_point2%point=face_zmax
      
      CALL get_cell_point_coordinate(cell_point1,point1)
      CALL get_cell_point_coordinate(cell_point2,point2)

! check whether the point defined by point1,point2 intersects the current triangle

      CALL line_triangle_intersection(point1,point2,triangle,	&
	                              intersection_point,intersection )
	       
      if (intersection) then
	
        CALL point_to_cell_no_checks(intersection_point,intersection_cell)
		
        mesh_face%cell%i=intersection_cell%i
        mesh_face%cell%j=intersection_cell%j
        mesh_face%cell%k=intersection_cell%k

! we must now find the closest normal to z face to the intersection point	
        cell_point1%cell%i=intersection_cell%i
        cell_point1%cell%j=intersection_cell%j
        cell_point1%cell%k=intersection_cell%k
        cell_point1%point=face_zmin
        CALL get_cell_point_coordinate(cell_point1,point1)
      
        cell_point2%cell%i=intersection_cell%i
        cell_point2%cell%j=intersection_cell%j
        cell_point2%cell%k=intersection_cell%k
        cell_point2%point=face_zmax
        CALL get_cell_point_coordinate(cell_point2,point2)
      
        dmin=xyz_distance(point1,intersection_point)
        dmax=xyz_distance(point2,intersection_point)
	
	if (dmin.lt.dmax) then
          mesh_face%point=face_zmin
	else
          mesh_face%point=face_zmax	
	end if
	
! calculate triangle normal and swap to the other side of the face if required to keep the mesh normal the same
        CALL triangle_normal(triangle,normal)
	if ( (normal%z.gt.0d0).AND.(mesh_face%point.EQ.face_zmin) ) then
          mesh_face%cell%k=mesh_face%cell%k-1
	  mesh_face%point=face_zmax
	else if ( (normal%z.lt.0d0).AND.(mesh_face%point.EQ.face_zmax) ) then
          mesh_face%cell%k=mesh_face%cell%k+1
	  mesh_face%point=face_zmin
	end if
	
	CALL cell_point_inside_mesh(mesh_face,inside)
	
	if (inside) then ! set face in the surface mesh array
	  surface_mesh(mesh_face%cell%i,mesh_face%cell%j,mesh_face%cell%k,mesh_face%point)=1
	end if
	
      end if
      
    end do ! next z cell
  end do ! next y cell 

  RETURN
  
  END SUBROUTINE mesh_triangle
