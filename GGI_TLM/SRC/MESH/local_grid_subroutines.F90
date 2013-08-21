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
!SUBROUTINE local_grid_line
!SUBROUTINE local_grid_triangle
!SUBROUTINE local_grid_tet
!SUBROUTINE allocate_local_grid
!SUBROUTINE get_closest_local_grid_corner
!SUBROUTINE grid_corner_coordinate
!SUBROUTINE dist_local_grid_point_to_line
!SUBROUTINE dist_local_grid_point_to_triangle
!SUBROUTINE get_face_list_normal_to_x
!SUBROUTINE get_face_list_normal_to_y
!SUBROUTINE get_face_list_normal_to_z
!
! NAME
!     SUBROUTINE local_grid_line
!
! DESCRIPTION
!     local_grid_line:
!
!     mesh a line in the local grid.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE local_grid_line(x1_in,y1_in,z1_in,x2_in,y2_in,z2_in)
  
USE local_mesh

IMPLICIT NONE

  real*8		:: x1_in,y1_in,z1_in
  real*8		:: x2_in,y2_in,z2_in

! local variables
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  
  integer 		:: ix1,iy1,iz1
  integer 		:: ix2,iy2,iz2
  
  integer 		:: offset(3)
  
  integer		:: ix,iy,iz
  
  integer		:: n_local_points,local_point
  integer		:: n_local_line_segments
  
  integer,allocatable	:: local_point_list(:,:)
  
  integer		:: point
  
! START

!  write(*,*)'CALLED local_grid_line'

! work out the ordering of the points i.e. which way to traverse the line 
! Ordering the points should make the mesh generation more consistent

  CALL order_points_2(x1_in,y1_in,z1_in,x2_in,y2_in,z2_in,x1,y1,z1,x2,y2,z2)

! get the line end points

  CALL get_closest_local_grid_corner(x1,y1,z1,ix1,iy1,iz1)
  CALL get_closest_local_grid_corner(x2,y2,z2,ix2,iy2,iz2)
  
!  write(*,*)'End point 1',ix1,iy1,iz1
!  write(*,*)'End point 2',ix2,iy2,iz2
  
! get the stepping direction in the x, y and z directions, ox,oy and oz  

  if (ix2.lt.ix1) then
    offset(1)=-1
  else if (ix2.gt.ix1) then
    offset(1)=+1
  else
    offset(1)=0
  end if

  if (iy2.lt.iy1) then
    offset(2)=-1
  else if (iy2.gt.iy1) then
    offset(2)=+1
  else
    offset(2)=0
  end if

  if (iz2.lt.iz1) then
    offset(3)=-1
  else if (iz2.gt.iz1) then
    offset(3)=+1
  else
    offset(3)=0
  end if

! set end points
  local_grid(ix1,iy1,iz1,local_corner)=1
  local_grid(ix2,iy2,iz2,local_corner)=1
  
! set the line between the end points in a crude way before iteratively improving the fit to the edge
  n_local_line_segments=abs(ix2-ix1)+abs(iy2-iy1)+abs(iz2-iz1)
  n_local_points=n_local_line_segments+1
  
  if (n_local_points.EQ.2) then
  
! two points and one edge only  
    CALL add_edge(ix1,iy1,iz1,ix2,iy2,iz2)
    
  else if (n_local_points.GT.2) then
  
    ALLOCATE( local_point_list(1:n_local_points,1:3) )

    local_point=1
    local_point_list(local_point,1)=ix1
    local_point_list(local_point,2)=iy1
    local_point_list(local_point,3)=iz1

! loop along line in x, y then z, setting points
    if (offset(1).NE.0) then
      iy=iy1
      iz=iz1
      do ix=ix1+offset(1),ix2,offset(1)
	local_point=local_point+1
	local_point_list(local_point,1)=ix
	local_point_list(local_point,2)=iy
	local_point_list(local_point,3)=iz
      end do
    end if
  
    if (offset(2).NE.0) then
      ix=ix2
      iz=iz1
      do iy=iy1+offset(2),iy2,offset(2)
	local_point=local_point+1
	local_point_list(local_point,1)=ix
	local_point_list(local_point,2)=iy
	local_point_list(local_point,3)=iz
      end do
    end if
  
    if (offset(3).NE.0) then
      ix=ix2
      iy=iy2
      do iz=iz1+offset(3),iz2,offset(3)
	local_point=local_point+1
	local_point_list(local_point,1)=ix
	local_point_list(local_point,2)=iy
	local_point_list(local_point,3)=iz
      end do
    end if
    
    if (local_point.NE.n_local_points) then
      write(*,*)'ERROR in local_grid_line'
      write(*,*)'Internal point counting error'
      write(*,*)'local_point   =',local_point
      write(*,*)'n_local_points=',n_local_points
      STOP
    end if

! iteratively improve the fit to the edge

    CALL iterative_improvement_points(x1,y1,z1,x2,y2,z2,local_point_list,n_local_points)

! put the edges into the local mesh      
    do point=1,n_local_points

      ix=local_point_list(point,1)
      iy=local_point_list(point,2)
      iz=local_point_list(point,3)
      local_grid(ix,iy,iz,local_corner)=1
      
      if (point.GE.2) then
        ix1=local_point_list(point-1,1)
        iy1=local_point_list(point-1,2)
        iz1=local_point_list(point-1,3)
        CALL add_edge(ix1,iy1,iz1,ix,iy,iz)
      end if
      
    end do
    
    DEALLOCATE( local_point_list )

  end if ! n_local_points.GT.0   
   
  RETURN
    
  END SUBROUTINE local_grid_line
!
! NAME
!     SUBROUTINE local_grid_triangle
!
! DESCRIPTION
!     local_grid_triangle:
!
!     mesh a triangle in the local grid.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE local_grid_triangle(x1_in,y1_in,z1_in,x2_in,y2_in,z2_in,x3_in,y3_in,z3_in)
  
USE local_mesh

IMPLICIT NONE

  real*8		:: x1_in,y1_in,z1_in
  real*8		:: x2_in,y2_in,z2_in
  real*8		:: x3_in,y3_in,z3_in

! local variables

  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
  
  integer	:: ix,iy,iz
  integer 	:: iface
  integer	:: last_n_faces_set
  integer	:: iteration
  
  logical	:: error_found
    
! START

!  write(*,*)'CALLED local_grid_triangle'
!  write(*,*)x1,y1,z1
!  write(*,*)x2,y2,z2
!  write(*,*)x3,y3,z3

  CALL order_points_3(x1_in,y1_in,z1_in,x2_in,y2_in,z2_in,x3_in,y3_in,z3_in,	&
                      x1,y1,z1,x2,y2,z2,x3,y3,z3)
		      
! set the edges of the triangle in the grid  
  CALL local_grid_line(x1,y1,z1,x2,y2,z2)
  CALL local_grid_line(x2,y2,z2,x3,y3,z3)
  CALL local_grid_line(x3,y3,z3,x1,y1,z1)
  
! work out where the faces normal to x, y and z are i.e. idetify the column which holds the face
  CALL get_face_list_normal_to_x()
  CALL get_face_list_normal_to_y()
  CALL get_face_list_normal_to_z()
  
  tot_n_faces=n_x_faces+n_y_faces+n_z_faces 
  n_faces_set=0
  last_n_faces_set=0
  iteration=0

10  CONTINUE  

    iteration=iteration+1

    CALL set_n_edge_faces(3)
    CALL set_n_edge_faces(2)
    
    if ( (n_faces_set.NE.tot_n_faces).AND.(n_faces_set.EQ.last_n_faces_set) ) then
! the iteration has not set any more faces so we have a problem with the mesh generation here

      write(*,*)'Error in local_grid_triangle'
      write(*,*)'Unable to complete filling the triangle surfaces'
      write(*,*)'Total number of faces to set       :',tot_n_faces
      write(*,*)'Number of faces which have been set:',n_faces_set    
      
      write(*,*)'point 1:',x1,y1,z1
      write(*,*)'point 2:',x2,y2,z2
      write(*,*)'point 3:',x3,y3,z3
      
      do ix=local_ixmin,local_ixmax
        do iy=local_iymin,local_iymax
          do iz=local_izmin,local_izmax
      
            if (local_grid(ix,iy,iz,local_corner).NE.0) then
	      write(*,*)'Set point',ix,iy,iz
            end if
	
          end do
        end do
      end do
	 
      STOP
      
    end if
    
    last_n_faces_set=n_faces_set
    
    if (n_faces_set.NE.tot_n_faces) GOTO 10 ! another iteration of the fill process  
    
    CALL check_edges(error_found)   
    
    if (error_found) then
    
      write(*,*)'Error found in check_edges - there is at least one free edge in a triangle surface mesh'   
      STOP

    end if

!  write(*,*)'CALLING iterative_improvement'
  CALL iterative_improvement_faces(x1,y1,z1,x2,y2,z2,x3,y3,z3)
  
! deallocate all corners and edges in the local grid

  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
      
        local_grid(ix,iy,iz,local_corner)=0
        local_grid(ix,iy,iz,local_centre)=0
        local_grid(ix,iy,iz,local_xedge) =0
        local_grid(ix,iy,iz,local_yedge) =0
        local_grid(ix,iy,iz,local_zedge) =0
	
      end do
    end do
  end do

  if ( allocated(local_x_faces) ) deallocate(local_x_faces)
  if ( allocated(local_y_faces) ) deallocate(local_y_faces)
  if ( allocated(local_z_faces) ) deallocate(local_z_faces)
  
  RETURN
  
  END SUBROUTINE local_grid_triangle
!
! NAME
!     SUBROUTINE local_grid_tetrahedron
!
! DESCRIPTION
!     local_grid_tetrahedron:
!
!     mesh a tet in the local grid.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE local_grid_tetrahedron(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
  
USE local_mesh

IMPLICIT NONE

  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
  real*8		:: x4,y4,z4

! local variables

  integer ix,iy,iz,point
  integer iz_tet_min,iz_tet_max
   
! START
  
! set the four triangles constituting the tet outer surface

  CALL  local_grid_triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3)
  
! copy triangle to local_grid_tet      
  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
        do point=1,8
  	  if (local_grid(ix,iy,iz,point).NE.0) then
	    local_grid_tet(ix,iy,iz,point)=local_grid(ix,iy,iz,point)
	    local_grid(ix,iy,iz,point)=0
	  end if
        end do
      end do
    end do
  end do
      
  CALL  local_grid_triangle(x1,y1,z1,x2,y2,z2,x4,y4,z4)
  
! copy triangle to local_grid_tet        
  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
        do point=1,8
  	  if (local_grid(ix,iy,iz,point).NE.0) then
	    local_grid_tet(ix,iy,iz,point)=local_grid(ix,iy,iz,point)
	    local_grid(ix,iy,iz,point)=0
	  end if
        end do
      end do
    end do
  end do
  
  CALL  local_grid_triangle(x2,y2,z2,x3,y3,z3,x4,y4,z4)
  
! copy triangle to local_grid_tet        
  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
        do point=1,8
  	  if (local_grid(ix,iy,iz,point).NE.0) then
	    local_grid_tet(ix,iy,iz,point)=local_grid(ix,iy,iz,point)
	    local_grid(ix,iy,iz,point)=0
	  end if
        end do
      end do
    end do
  end do
  
  CALL  local_grid_triangle(x3,y3,z3,x1,y1,z1,x4,y4,z4)
  
! copy triangle to local_grid_tet        
  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
        do point=1,8
  	  if (local_grid(ix,iy,iz,point).NE.0) then
	    local_grid_tet(ix,iy,iz,point)=local_grid(ix,iy,iz,point)
	    local_grid(ix,iy,iz,point)=0
	  end if
        end do
      end do
    end do
  end do
  
! fill the cells by looping over the xy plane and filling columns of z cells

! loop over x and y
  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
    
      iz_tet_min=not_set
      iz_tet_max=not_set
    
! STAGE 1. look for a zface to be set at this ix,iy
      do iz=local_izmin,local_izmax
      
        if ( (local_grid_tet(ix,iy,iz,local_zface).NE.0).AND.(iz_tet_min.EQ.not_set) ) then  
! the zmin face of this cell is set
          iz_tet_min=iz
        end if 
	
      end do
      
! STAGE 2. look for a second zface to be set at this ix,iy
      if (iz_tet_min.NE.not_set) then
! we have found a face normal to z so look for a second face normal to z which bounds the tet cells
        do iz=iz_tet_min+1,local_izmax
      
          if ( (local_grid_tet(ix,iy,iz,local_zface).NE.0).AND.(iz_tet_max.EQ.not_set) ) then  
! the zmin face of this cell is set
            iz_tet_max=iz
          end if 
	
        end do ! next iz
	
      end if ! iz_tet_min.NE.not_set
     
! STAGE 3. set the cells between the found faces
      if ( (iz_tet_min.NE.not_set).AND.(iz_tet_max.NE.not_set) ) then
! we have found two faces which bound the tet cells so set all the cells in between the faces found
        do iz=iz_tet_min,iz_tet_max-1
      
          local_grid_tet(ix,iy,iz,local_centre)=1

        end do ! next iz
	
      end if ! two bounding faces found
      
    end do ! next loop over iy
  end do ! next loop over ix
  
  RETURN
  
  END SUBROUTINE local_grid_tetrahedron
!
! NAME
!     SUBROUTINE allocate_local_grid
!
! DESCRIPTION
!     allocate_local_grid:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE allocate_local_grid(min_point,max_point)
  
USE geometry_types
USE local_mesh

IMPLICIT NONE

  type(xyz)		:: min_point,max_point

! local variables

  type(ijk)		:: min_cell,max_cell
   
! START

  CALL point_to_cell_no_checks(min_point,min_cell)
  CALL point_to_cell_no_checks(max_point,max_cell)
  
! Define the dimensions of the local grid add an extra cell all the way round
  local_ixmin=min_cell%i-1
  local_iymin=min_cell%j-1
  local_izmin=min_cell%k-1
  
  local_ixmax=max_cell%i+1
  local_iymax=max_cell%j+1
  local_izmax=max_cell%k+1
  
  local_xmin=mesh_xmin+(local_ixmin-1)*dl
  local_ymin=mesh_ymin+(local_iymin-1)*dl
  local_zmin=mesh_zmin+(local_izmin-1)*dl
    
  ALLOCATE( local_grid(local_ixmin:local_ixmax,local_iymin:local_iymax,local_izmin:local_izmax,1:8) )
  
  local_grid(local_ixmin:local_ixmax,local_iymin:local_iymax,local_izmin:local_izmax,1:8)=0

  RETURN
  
  END SUBROUTINE allocate_local_grid
!
! NAME
!     SUBROUTINE allocate_local_grid_tet
!
! DESCRIPTION
!     allocate_local_grid_tet:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE allocate_local_grid_tet(min_point,max_point)
  
USE geometry_types
USE local_mesh

IMPLICIT NONE

  type(xyz)		:: min_point,max_point

! local variables

  type(ijk)		:: min_cell,max_cell
   
! START

  CALL point_to_cell_no_checks(min_point,min_cell)
  CALL point_to_cell_no_checks(max_point,max_cell)
  
! Define the dimensions of the local grid add an extra cell all the way round
  local_ixmin=min_cell%i-1
  local_iymin=min_cell%j-1
  local_izmin=min_cell%k-1
  
  local_ixmax=max_cell%i+1
  local_iymax=max_cell%j+1
  local_izmax=max_cell%k+1
  
  local_xmin=mesh_xmin+(local_ixmin-1)*dl
  local_ymin=mesh_ymin+(local_iymin-1)*dl
  local_zmin=mesh_zmin+(local_izmin-1)*dl
      
  ALLOCATE( local_grid_tet(local_ixmin:local_ixmax,local_iymin:local_iymax,local_izmin:local_izmax,1:8) )
  
  local_grid_tet(local_ixmin:local_ixmax,local_iymin:local_iymax,local_izmin:local_izmax,1:8)=0

  RETURN
  
  END SUBROUTINE allocate_local_grid_tet
!
! NAME
!     SUBROUTINE get_closest_local_grid_corner
!
! DESCRIPTION
!     get_closest_local_grid_corner:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE get_closest_local_grid_corner(x,y,z,ix,iy,iz)
  
USE local_mesh

IMPLICIT NONE

  real*8	:: x,y,z
  integer	:: ix,iy,iz

! local variables
   
! START

  ix=nint((x-mesh_xmin)/dl)+1
  iy=nint((y-mesh_ymin)/dl)+1
  iz=nint((z-mesh_zmin)/dl)+1
  
  if ((ix.lt.local_ixmin).OR.(ix.gt.local_ixmax)) goto 9000
  if ((iy.lt.local_iymin).OR.(iy.gt.local_iymax)) goto 9000
  if ((iz.lt.local_izmin).OR.(iz.gt.local_izmax)) goto 9000

  RETURN

9000 write(*,*)'Error in get_closest_local_grid_corner'
     write(*,*)'The cell is outside the defined local grid'
     write(*,*)'x =',x,' y =',y,' z =',z
     write(*,*)'ix =',ix,' iy =',iy,' iz =',iz
     write(*,*)'ixmin =',local_ixmin ,'iymin =',local_iymin ,'izmin =',local_izmin
     write(*,*)'ixmax =',local_ixmax ,'iymax =',local_iymax ,'izmax =',local_izmax
     STOP
  
  END SUBROUTINE get_closest_local_grid_corner
!
! NAME
!     SUBROUTINE grid_corner_coordinate
!
! DESCRIPTION
!     grid_corner_coordinate
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE grid_corner_coordinate(ix,iy,iz,x,y,z)
  
USE local_mesh

IMPLICIT NONE

  integer	:: ix,iy,iz
  real*8	:: x,y,z

! local variables


! START

  x=mesh_xmin+(ix-1)*dl
  y=mesh_ymin+(iy-1)*dl
  z=mesh_zmin+(iz-1)*dl
  
  RETURN
  
  END SUBROUTINE grid_corner_coordinate
!
! NAME
!     SUBROUTINE dist_local_grid_point_to_line
!
! DESCRIPTION
!     dist_local_grid_point_to_line:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE dist_local_grid_point_to_line(ix,iy,iz,x1,y1,z1,x2,y2,z2,dist)
  
USE local_mesh

IMPLICIT NONE

  integer	:: ix,iy,iz
  real*8	:: x1,y1,z1,x2,y2,z2,dist

! local variables

  real*8	:: x,y,z
  
  real*8	:: vx,vy,vz,length
  real*8	:: t
  real*8	:: px,py,pz
  real*8	:: ux,uy,uz
   
! START

! get the coordinates of the logal grid point

  CALL grid_corner_coordinate(ix,iy,iz,x,y,z)

! get the vector along the line
  vx=x2-x1
  vy=y2-y1
  vz=z2-z1
  
  length=sqrt(vx*vx+vy*vy+vz*vz)
  
  if (length.eq.0d0) then
! line length 0 so return distance of point to end point 1
    dist=sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)
    RETURN
  end if
  
! the equation of the line is p=end1+v*t

! get the vector from end point 1 to the test point
  px=x-x1
  py=y-y1
  pz=z-z1
  
! get t the normalised distance along the (inifinite) line of the closest point

  t=(px*vx+py*vy+pz*vz)/(vx*vx+vy*vy+vz*vz)  
  
  if (t.lt.0d0) then
! point is closest to end point 1 of the finite line segment

    dist=sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)
    RETURN
    
  else if (t.gt.1d0) then
! point is closest to end point 2 of the finite line segment

    dist=sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2)
    RETURN
    
  else
  
    ux=x1+t*vx
    uy=y1+t*vy
    uz=z1+t*vz
    dist=sqrt((x-ux)**2+(y-uy)**2+(z-uz)**2)
  
  end if

  RETURN
  
  END SUBROUTINE dist_local_grid_point_to_line
!
! NAME
!     SUBROUTINE dist_local_grid_point_to_triangle
!
! DESCRIPTION
!     dist_local_grid_point_to_triangle:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/5/2013 CJS 
!
SUBROUTINE dist_local_grid_point_to_triangle(ix,iy,iz,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist)
  
USE local_mesh

IMPLICIT NONE

  integer	:: ix,iy,iz
  real*8	:: x1,y1,z1,x2,y2,z2,x3,y3,z3,dist

! local variables

  real*8	:: x,y,z
  
  real*8	:: v1x,v1y,v1z
  real*8	:: v2x,v2y,v2z
  real*8	:: length
  real*8	:: normx,normy,normz

  real*8	:: px,py,pz
  
  real*8	:: tx,ty,tz
  
  real*8	:: a,b,c,d,e,f
  
  real*8	:: s,t
   
! START

! get the coordinates of the logal grid point

  CALL grid_corner_coordinate(ix,iy,iz,x,y,z)

! get two vectors in the plane of the triangle
  v1x=x2-x1
  v1y=y2-y1
  v1z=z2-z1
  
  v2x=x3-x1
  v2y=y3-y1
  v2z=z3-z1

! get the normal vector n=v1^v2

  normx=v1y*v2z-v1z*v2y
  normy=v1z*v2x-v1x*v2z
  normz=v1x*v2y-v1y*v2x

  length=sqrt(normx*normx+normy*normy+normz*normz)
  
  if (length.eq.0d0) then
! triangle area =0 so return distance of point to line from point 1 to point 2
    CALL dist_local_grid_point_to_line(ix,iy,iz,x1,y1,z1,x2,y2,z2,dist)
    RETURN
  end if
  
  normx=normx/length
  normy=normy/length
  normz=normz/length
  
! the equation of the line is R.n=p

! get the vector from point 1 to the test point
  px=x-x1
  py=y-y1
  pz=z-z1

! distance from plane to point is p.n  
  dist=abs(px*normx+py*normy+pz*normz)
  
! work out the projection of the point onto the plane of the triangle

  tx=px-dist*normx
  ty=py-dist*normy
  tz=pz-dist*normz
    
! the point t in the triangle is expressed as t=s(v1)+t(v2)

  a=v1x*v1x+v1y*v1y+v1z*v1z
  b=v1x*v2x+v1y*v2y+v1z*v2z
  c=v2x*v2x+v2y*v2y+v2z*v2z
  d=-(v1x*tx+v1y*ty+v1z*tz)
  e=-(v2x*tx+v2y*ty+v2z*tz)
  f= (tx*tx+ty*ty+tz*tz)
  
! work out the parametric coordinates

  s=(b*e-c*d)/(a*c-b*b)
  t=(b*d-a*e)/(a*c-b*b)
  
  if ( (s.GE.0d0).and.(t.GE.0d0).AND.(s+t.LT.1d0) ) then
! inside triangle so return the distance of the point from the plane
    RETURN
  else if ( (s+t.GT.1d0) ) then
    CALL dist_local_grid_point_to_line(ix,iy,iz,x2,y2,z2,x3,y3,z3,dist)
    RETURN
  else if ( (s.LT.0d0) ) then
    CALL dist_local_grid_point_to_line(ix,iy,iz,x1,y1,z1,x3,y3,z3,dist)
    RETURN
  else if ( (t.LT.0d0) ) then
    CALL dist_local_grid_point_to_line(ix,iy,iz,x1,y1,z1,x2,y2,z2,dist)
    RETURN
  end if

  RETURN
  
  END SUBROUTINE dist_local_grid_point_to_triangle
!
! NAME
!     SUBROUTINE get_face_list_normal_to_x
!
! DESCRIPTION
!     get_face_list_normal_to_x:
!
!     work out the number of triangle faces normal to x and identify the y and z columns
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/5/2013 CJS 
!
SUBROUTINE get_face_list_normal_to_x()
  
USE local_mesh

IMPLICIT NONE

! local variables

  integer loop
  integer ix,iy,iz
  integer iz_tri_min,iz_tri_max
  
  logical found_edge
  logical found_first
  logical found_second
   
! START

! two loops, first to count the faces and second to fill the face list
  
  do loop=1,2

    n_x_faces=0
    
! loop over the yz plane and filling columns of x cells

! loop over y
    do iy=local_iymin,local_iymax

! loop over z looking for triangle edges in the y direction

      found_first =.FALSE.
      found_second=.FALSE.
      
      do iz=local_izmin,local_izmax
    
! loop over x looking for an y directed edge

        found_edge=.FALSE.

        do ix=local_ixmin,local_ixmax
	
	  if (local_grid(ix,iy,iz,local_yedge).NE.0) then
	    found_edge=.TRUE.
	  end if
	
	end do ! next ix
	
	if (found_edge) then
	
	  if (.NOT.found_first) then
	    found_first =.TRUE.
	    iz_tri_min=iz
	  else if (.NOT.found_second) then
	    found_second=.TRUE.
	    iz_tri_max=iz-1
	  else
	    write(*,*)'Error in get_face_list_normal_to_x'
	    write(*,*)'Found more than 2 edges, iz=',iz
	    STOP
	  end if
	  
	end if ! found_edge

      end do ! next iz
      
      if (found_first.AND.found_second) then
	
	do iz=iz_tri_min,iz_tri_max
	
	  n_x_faces=n_x_faces+1
	
	  if (loop.eq.2) then
	  
	    local_x_faces(n_x_faces,1)=not_set
	    local_x_faces(n_x_faces,2)=iy
	    local_x_faces(n_x_faces,3)=iz
	    
	  end if
	  
	end do ! next iz inside triangle
	
      end if
      
    end do ! next iy
    	
    if (loop.eq.1) then

      if (n_x_faces.gt.0) then
        ALLOCATE ( local_x_faces(1:n_x_faces,1:3) )
        local_x_faces(1:n_x_faces,1:3)=0
      end if
      
    end if

  end do ! next loop

  
  RETURN
  
  END SUBROUTINE get_face_list_normal_to_x
!
! NAME
!     SUBROUTINE get_face_list_normal_to_y
!
! DESCRIPTION
!     get_face_list_normal_to_y:
!
!     work out the number of triangle faces normal to y and identify the z and x columns
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/5/2013 CJS 
!
SUBROUTINE get_face_list_normal_to_y()
  
USE local_mesh

IMPLICIT NONE

! local variables

  integer loop
  integer ix,iy,iz
  integer ix_tri_min,ix_tri_max
  
  logical found_edge
  logical found_first
  logical found_second
   
! START

! two loops, first to count the faces and second to fill the face list
  
  do loop=1,2

    n_y_faces=0
    
! loop over the zx plane and filling columns of y cells

! loop over z
    do iz=local_izmin,local_izmax

! loop over x looking for triangle edges in the z direction

      found_first =.FALSE.
      found_second=.FALSE.
      
      do ix=local_ixmin,local_ixmax
    
! loop over y looking for a z directed edge

        found_edge=.FALSE.

        do iy=local_iymin,local_iymax
	
	  if (local_grid(ix,iy,iz,local_zedge).NE.0) then
	    found_edge=.TRUE.
	  end if
	
	end do ! next iy
	
	if (found_edge) then
	
	  if (.NOT.found_first) then
	    found_first =.TRUE.
	    ix_tri_min=ix
	  else if (.NOT.found_second) then
	    found_second=.TRUE.
	    ix_tri_max=ix-1
	  else
	    write(*,*)'Error in get_face_list_normal_to_y'
	    write(*,*)'Found more than 2 edges, ix=',ix
	    STOP
	  end if
	  
	end if ! found_edge

      end do ! next ix
      
      if (found_first.AND.found_second) then
	
	do ix=ix_tri_min,ix_tri_max
	
	  n_y_faces=n_y_faces+1
	
	  if (loop.eq.2) then
	  
	    local_y_faces(n_y_faces,1)=ix
	    local_y_faces(n_y_faces,2)=not_set
	    local_y_faces(n_y_faces,3)=iz
	    
	  end if
	  
	end do ! next ix inside triangle
	
      end if
      
    end do ! next iz
    	
    if (loop.eq.1) then

      if (n_y_faces.gt.0) then
        ALLOCATE ( local_y_faces(1:n_y_faces,1:3) )
        local_y_faces(1:n_y_faces,1:3)=0
      end if
      
    end if

  end do ! next loop

  
  RETURN
  
  END SUBROUTINE get_face_list_normal_to_y

!
! NAME
!     SUBROUTINE get_face_list_normal_to_z
!
! DESCRIPTION
!     get_face_list_normal_to_z:
!
!     work out the number of triangle faces normal to z and identify the x and y columns
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/5/2013 CJS 
!
SUBROUTINE get_face_list_normal_to_z()
  
USE local_mesh

IMPLICIT NONE

! local variables

  integer loop
  integer ix,iy,iz
  integer iy_tri_min,iy_tri_max
  
  logical found_edge
  logical found_first
  logical found_second
   
! START

! two loops, first to count the faces and second to fill the face list
  
  do loop=1,2

    n_z_faces=0
    
! loop over the xy plane and filling columns of z cells

! loop over x 
    do ix=local_ixmin,local_ixmax

! loop over y looking for triangle edges in the x direction

      found_first =.FALSE.
      found_second=.FALSE.
      
      do iy=local_iymin,local_iymax
    
! loop over z looking for an x directed edge

        found_edge=.FALSE.

        do iz=local_izmin,local_izmax
	
	  if (local_grid(ix,iy,iz,local_xedge).NE.0) then
	    found_edge=.TRUE.
	  end if
	
	end do ! next iz
	
	if (found_edge) then
	
	  if (.NOT.found_first) then
	    found_first =.TRUE.
	    iy_tri_min=iy
	  else if (.NOT.found_second) then
	    found_second=.TRUE.
	    iy_tri_max=iy-1
	  else
	    write(*,*)'Error in get_face_list_normal_to_z'
	    write(*,*)'Found more than 2 edges, iy=',iy
	    STOP
	  end if
	  
	end if ! found_edge

      end do ! next iy
      
      if (found_first.AND.found_second) then
	
	do iy=iy_tri_min,iy_tri_max
	
	  n_z_faces=n_z_faces+1
	
	  if (loop.eq.2) then
	  
	    local_z_faces(n_z_faces,1)=ix
	    local_z_faces(n_z_faces,2)=iy
	    local_z_faces(n_z_faces,3)=not_set
	    
	  end if
	  
	end do ! next iy inside triangle
	
      end if
      
    end do ! next ix
    	
    if (loop.eq.1) then

      if (n_z_faces.gt.0) then
        ALLOCATE ( local_z_faces(1:n_z_faces,1:3) )
        local_z_faces(1:n_z_faces,1:3)=0
      end if
      
    end if

  end do ! next loop

  
  RETURN
  
  END SUBROUTINE get_face_list_normal_to_z
!
! NAME
!     SUBROUTINE subtract_face
!
! DESCRIPTION
!     remove a face and associated edges from the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 21/5/2013 CJS 
!
SUBROUTINE subtract_face(ix,iy,iz,face)
  
USE local_mesh

IMPLICIT NONE

  integer ix,iy,iz,face

! local variables
   
! START

! subtract the face

!  write(*,*)'Subtract face',ix,iy,iz,face

  if (face.eq.local_xface) then
  
    local_grid(ix  ,iy  ,iz  ,local_xface)=0
		
! correct the associated edge data
    local_grid(ix  ,iy  ,iz  ,local_yedge)=local_grid(ix  ,iy  ,iz  ,local_yedge)-1
    local_grid(ix  ,iy  ,iz+1,local_yedge)=local_grid(ix  ,iy  ,iz+1,local_yedge)-1
    local_grid(ix  ,iy  ,iz  ,local_zedge)=local_grid(ix  ,iy  ,iz  ,local_zedge)-1
    local_grid(ix  ,iy+1,iz  ,local_zedge)=local_grid(ix  ,iy+1,iz  ,local_zedge)-1

  else if (face.eq.local_yface) then
  
    local_grid(ix  ,iy  ,iz  ,local_yface)=0
		
! correct the associated edge data
    local_grid(ix  ,iy  ,iz  ,local_xedge)=local_grid(ix  ,iy  ,iz  ,local_xedge)-1
    local_grid(ix  ,iy  ,iz+1,local_xedge)=local_grid(ix  ,iy  ,iz+1,local_xedge)-1
    local_grid(ix  ,iy  ,iz  ,local_zedge)=local_grid(ix  ,iy  ,iz  ,local_zedge)-1
    local_grid(ix+1,iy  ,iz  ,local_zedge)=local_grid(ix+1,iy  ,iz  ,local_zedge)-1

  else if (face.eq.local_zface) then
  
    local_grid(ix  ,iy  ,iz  ,local_zface)=0
		
! correct the associated edge data
    local_grid(ix  ,iy  ,iz  ,local_xedge)=local_grid(ix  ,iy  ,iz  ,local_xedge)-1
    local_grid(ix  ,iy+1,iz  ,local_xedge)=local_grid(ix  ,iy+1,iz  ,local_xedge)-1
    local_grid(ix  ,iy  ,iz  ,local_yedge)=local_grid(ix  ,iy  ,iz  ,local_yedge)-1
    local_grid(ix+1,iy  ,iz  ,local_yedge)=local_grid(ix+1,iy  ,iz  ,local_yedge)-1
    
  end if
  
END SUBROUTINE subtract_face
!
! NAME
!     SUBROUTINE add_face
!
! DESCRIPTION
!     add a face and associated edges to the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 21/5/2013 CJS 
!
SUBROUTINE add_face(ix,iy,iz,face)
  
USE local_mesh

IMPLICIT NONE

  integer ix,iy,iz,face

! local variables
   
! START

! add the face
!  write(*,*)'Add face',ix,iy,iz,face

  if (face.eq.local_xface) then
  
    local_grid(ix  ,iy  ,iz  ,local_xface)=1
		
! increment the associated edge data
    local_grid(ix  ,iy  ,iz  ,local_yedge)=local_grid(ix  ,iy  ,iz  ,local_yedge)+1
    local_grid(ix  ,iy  ,iz+1,local_yedge)=local_grid(ix  ,iy  ,iz+1,local_yedge)+1
    local_grid(ix  ,iy  ,iz  ,local_zedge)=local_grid(ix  ,iy  ,iz  ,local_zedge)+1
    local_grid(ix  ,iy+1,iz  ,local_zedge)=local_grid(ix  ,iy+1,iz  ,local_zedge)+1
		
! make sure that the point data is set
    local_grid(ix  ,iy  ,iz  ,local_corner)=1
    local_grid(ix  ,iy+1,iz  ,local_corner)=1
    local_grid(ix  ,iy+1,iz+1,local_corner)=1
    local_grid(ix  ,iy  ,iz+1,local_corner)=1

  else if (face.eq.local_yface) then
  
    local_grid(ix  ,iy  ,iz  ,local_yface)=1
		
! increment the associated edge data
    local_grid(ix  ,iy  ,iz  ,local_xedge)=local_grid(ix  ,iy  ,iz  ,local_xedge)+1
    local_grid(ix  ,iy  ,iz+1,local_xedge)=local_grid(ix  ,iy  ,iz+1,local_xedge)+1
    local_grid(ix  ,iy  ,iz  ,local_zedge)=local_grid(ix  ,iy  ,iz  ,local_zedge)+1
    local_grid(ix+1,iy  ,iz  ,local_zedge)=local_grid(ix+1,iy  ,iz  ,local_zedge)+1
		
! make sure that the point data is set
    local_grid(ix  ,iy  ,iz  ,local_corner)=1
    local_grid(ix+1,iy  ,iz  ,local_corner)=1
    local_grid(ix+1,iy  ,iz+1,local_corner)=1
    local_grid(ix  ,iy  ,iz+1,local_corner)=1

  else if (face.eq.local_zface) then
  
    local_grid(ix  ,iy  ,iz  ,local_zface)=1
		
! increment the associated edge data
    local_grid(ix  ,iy  ,iz  ,local_xedge)=local_grid(ix  ,iy  ,iz  ,local_xedge)+1
    local_grid(ix  ,iy+1,iz  ,local_xedge)=local_grid(ix  ,iy+1,iz  ,local_xedge)+1
    local_grid(ix  ,iy  ,iz  ,local_yedge)=local_grid(ix  ,iy  ,iz  ,local_yedge)+1
    local_grid(ix+1,iy  ,iz  ,local_yedge)=local_grid(ix+1,iy  ,iz  ,local_yedge)+1
		
! make sure that the point data is set
    local_grid(ix  ,iy  ,iz  ,local_corner)=1
    local_grid(ix  ,iy+1,iz  ,local_corner)=1
    local_grid(ix+1,iy+1,iz  ,local_corner)=1
    local_grid(ix+1,iy  ,iz  ,local_corner)=1
    
  end if
  
END SUBROUTINE add_face
!
! NAME
!     SUBROUTINE get_n_edges
!
! DESCRIPTION
!     add a face and associated edges to the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 21/5/2013 CJS 
!
SUBROUTINE get_n_edges(ix,iy,iz,face,value,n_edges)
  
USE local_mesh

IMPLICIT NONE

  integer ix,iy,iz,face,value,n_edges

! local variables
   
! START

  n_edges=0
  
  if      (face.eq.local_xface) then
  
    if (local_grid(ix  ,iy  ,iz  ,local_yedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix  ,iy  ,iz+1,local_yedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix  ,iy  ,iz  ,local_zedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix  ,iy+1,iz  ,local_zedge).EQ.value) n_edges=n_edges+1

  else if (face.eq.local_yface) then

    if (local_grid(ix  ,iy  ,iz  ,local_xedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix  ,iy  ,iz+1,local_xedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix  ,iy  ,iz  ,local_zedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix+1,iy  ,iz  ,local_zedge).EQ.value) n_edges=n_edges+1

  else if (face.eq.local_zface) then

    if (local_grid(ix  ,iy  ,iz  ,local_xedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix  ,iy+1,iz  ,local_xedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix  ,iy  ,iz  ,local_yedge).EQ.value) n_edges=n_edges+1
    if (local_grid(ix+1,iy  ,iz  ,local_yedge).EQ.value) n_edges=n_edges+1
    
  end if
  
  RETURN
  
END SUBROUTINE get_n_edges
!
! NAME
!     SUBROUTINE get_n_edges
!
! DESCRIPTION
!     add a face and associated edges to the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 21/5/2013 CJS 
!
SUBROUTINE get_n_faces(ix,iy,iz,n_faces)
  
USE local_mesh

IMPLICIT NONE

  integer ix,iy,iz,n_faces

! local variables
   
! START

  n_faces= local_grid(ix  ,iy  ,iz  ,local_xface)	&
          +local_grid(ix  ,iy  ,iz  ,local_yface)	&
          +local_grid(ix  ,iy  ,iz  ,local_zface)	&
          +local_grid(ix+1,iy  ,iz  ,local_xface)	&
          +local_grid(ix  ,iy+1,iz  ,local_yface)	&
          +local_grid(ix  ,iy  ,iz+1,local_zface)
  
  RETURN
  
END SUBROUTINE get_n_faces

!
! NAME
!     SUBROUTINE set_n_edge_faces
!
! DESCRIPTION
!
!     set_three_n_faces
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 28/5/2013 CJS 
!
SUBROUTINE set_n_edge_faces(n_edges_set)
  
USE local_mesh

IMPLICIT NONE

  integer n_edges_set

! local variables

  integer iface
  integer ix,iy,iz
  
  integer n_edges_value_1
  integer n_edges_value_2
     
! START

! reset all cell centres

! project across in all directions setting faces which have at least three edges set
  do iface=1,n_x_faces
  
    ix=local_x_faces(iface,1)
    iy=local_x_faces(iface,2)
    iz=local_x_faces(iface,3)
    
    if (ix.EQ.not_set) then
    
      do ix=local_ixmin,local_ixmax
      
        CALL get_n_edges(ix,iy,iz,local_xface,1,n_edges_value_1)
        CALL get_n_edges(ix,iy,iz,local_xface,2,n_edges_value_2)
	
	if ( (n_edges_value_1.GE.n_edges_set).AND.(n_edges_value_2.EQ.0).AND.(local_x_faces(iface,1).EQ.not_set) ) then
	
          CALL  add_face(ix,iy,iz,local_xface)
	  n_faces_set=n_faces_set+1
	  local_x_faces(iface,1)=ix
	  
	end if
	
      end do
      
    end if ! face not set
    
  end do
  
  do iface=1,n_y_faces
  
    ix=local_y_faces(iface,1)
    iy=local_y_faces(iface,2)
    iz=local_y_faces(iface,3)
    
    if (iy.EQ.not_set) then
    
      do iy=local_iymin,local_iymax
      
        CALL get_n_edges(ix,iy,iz,local_yface,1,n_edges_value_1)
        CALL get_n_edges(ix,iy,iz,local_yface,2,n_edges_value_2)
	
	if ( (n_edges_value_1.GE.n_edges_set).AND.(n_edges_value_2.EQ.0).AND.(local_y_faces(iface,2).EQ.not_set) ) then
	
          CALL  add_face(ix,iy,iz,local_yface)
	  n_faces_set=n_faces_set+1
	  local_y_faces(iface,2)=iy
	  
	end if

      end do
      
    end if ! not_set
    
  end do
  
  do iface=1,n_z_faces
  
    ix=local_z_faces(iface,1)
    iy=local_z_faces(iface,2)
    iz=local_z_faces(iface,3)
    
    if (iz.EQ.not_set) then
      
      do iz=local_izmin,local_izmax
      
        CALL get_n_edges(ix,iy,iz,local_zface,1,n_edges_value_1)
        CALL get_n_edges(ix,iy,iz,local_zface,2,n_edges_value_2)
	
	if ( (n_edges_value_1.GE.n_edges_set).AND.(n_edges_value_2.EQ.0).AND.(local_z_faces(iface,3).EQ.not_set) ) then
	
          CALL  add_face(ix,iy,iz,local_zface)
	  n_faces_set=n_faces_set+1
	  local_z_faces(iface,3)=iz
	  
	end if

      end do
      
    end if ! not_set
    
  end do

  
  RETURN
  
  END SUBROUTINE set_n_edge_faces
!
! NAME
!     SUBROUTINE iterative_improvement_faces
!
! DESCRIPTION
!
!     
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/5/2013 CJS 
!
SUBROUTINE iterative_improvement_faces(x1,y1,z1,x2,y2,z2,x3,y3,z3)
  
USE local_mesh

IMPLICIT NONE
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3

! local variables

  integer ix,iy,iz
  
  integer n_moved,n_faces
     
  real*8	:: dist_1,dist_2
  
  integer,parameter	:: max_iterations=1000000
  integer		:: n_iterations
  
     
! START
	    
  n_iterations=0

10  CONTINUE

    n_moved=0
    n_iterations=n_iterations+1
    
    if (n_iterations.GT.max_iterations) then
      write(*,*)'Maximum number of iterations exceeded in iterative_improvement_points'
      STOP
    end if
    
! loop over the triangle mesh looking for 'corner' cells in which 3 faces are set and 
! flip the faces if this improves the fit between the meshed surface and the original triangle
      
    do ix=local_ixmin,local_ixmax-1
      do iy=local_iymin,local_iymax-1
        do iz=local_izmin,local_izmax-1
	
          CALL get_n_faces(ix,iy,iz,n_faces)
	  
  	  if (n_faces.eq.3) then
    	    
! work out the point at which the faces intersect and test whether meshing the opposite faces provides
! an improved fit to the original triangle   

! check point ix  ,iy  ,iz
            if (       local_grid(ix  ,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix  ,iy  ,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix+1,iy+1,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)

!  point ix+1,iy  ,iz
          
	     else if ( local_grid(ix+1,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix+1,iy  ,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix  ,iy+1,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)
            
!  point ix+1,iy+1,iz
          
	     else if ( local_grid(ix+1,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy+1,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix+1,iy+1,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix  ,iy  ,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)
  
!  point ix  ,iy+1,iz
          
	     else if ( local_grid(ix  ,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy+1,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix  ,iy+1,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix+1,iy  ,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)

!  point ix  ,iy  ,iz+1
          
	     else if ( local_grid(ix  ,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz+1,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix  ,iy  ,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix+1,iy+1,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)
  
!  point ix+1,iy  ,iz+1
          
	     else if ( local_grid(ix+1,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy  ,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz+1,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix+1,iy  ,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix  ,iy+1,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)
  
!  point ix+1,iy+1,iz+1
          
	     else if ( local_grid(ix+1,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy+1,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz+1,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix+1,iy+1,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix  ,iy  ,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)
 
!  point ix  ,iy+1,iz+1	    
          
	     else if ( local_grid(ix  ,iy  ,iz  ,local_xface)	&
                      +local_grid(ix  ,iy+1,iz  ,local_yface)	&
                      +local_grid(ix  ,iy  ,iz+1,local_zface).EQ.3) then
	       CALL dist_local_grid_point_to_triangle(ix  ,iy+1,iz+1,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_1)
	       CALL dist_local_grid_point_to_triangle(ix+1,iy  ,iz  ,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist_2)
	    
	     else
	       write(*,*)'Error in iterative_improvement'
	       write(*,*)'Cannot find 3 face point'
!	       STOP

	     end if 
	     
	     if (dist_2.LT.dist_1) then
	       n_moved=n_moved+1
	       CALL swap_faces(ix,iy,iz)
	     end if
	      	    
  	  end if ! n_faces.eq.3

        end do
      end do
    end do

    if ((n_moved.NE.0)) GOTO 10
  
  RETURN
  
  END SUBROUTINE iterative_improvement_faces
!
! NAME
!     SUBROUTINE swap_faces
!
! DESCRIPTION
!     swap faces across the cell
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/5/2013 CJS 
!
SUBROUTINE swap_faces(ix,iy,iz)
  
USE local_mesh

IMPLICIT NONE

  integer ix,iy,iz,n_faces

! local variables
   
! START

! x faces
  if      ((local_grid(ix  ,iy  ,iz  ,local_xface).EQ.1).AND.(local_grid(ix+1,iy  ,iz  ,local_xface).EQ.0)) then
    CALL subtract_face(ix  ,iy  ,iz  ,local_xface)
    CALL      add_face(ix+1,iy  ,iz  ,local_xface)
  else if ((local_grid(ix  ,iy  ,iz  ,local_xface).EQ.0).AND.(local_grid(ix+1,iy  ,iz  ,local_xface).EQ.1)) then
    CALL subtract_face(ix+1,iy  ,iz  ,local_xface)
    CALL      add_face(ix  ,iy  ,iz  ,local_xface)  
  end if

! y faces
  if      ((local_grid(ix  ,iy  ,iz  ,local_yface).EQ.1).AND.(local_grid(ix  ,iy+1,iz  ,local_yface).EQ.0)) then
    CALL subtract_face(ix  ,iy  ,iz  ,local_yface)
    CALL      add_face(ix  ,iy+1,iz  ,local_yface) 
  else if ((local_grid(ix  ,iy  ,iz  ,local_yface).EQ.0).AND.(local_grid(ix  ,iy+1,iz  ,local_yface).EQ.1)) then
    CALL subtract_face(ix  ,iy+1,iz  ,local_yface)
    CALL      add_face(ix  ,iy  ,iz  ,local_yface)   
  end if

! z faces
  if      ((local_grid(ix  ,iy  ,iz  ,local_zface).EQ.1).AND.(local_grid(ix  ,iy  ,iz+1,local_zface).EQ.0)) then
    CALL subtract_face(ix  ,iy  ,iz  ,local_zface)
    CALL      add_face(ix  ,iy  ,iz+1,local_zface)   
  else if ((local_grid(ix  ,iy  ,iz  ,local_zface).EQ.0).AND.(local_grid(ix  ,iy  ,iz+1,local_zface).EQ.1)) then
    CALL subtract_face(ix  ,iy  ,iz+1,local_zface)
    CALL      add_face(ix  ,iy  ,iz  ,local_zface)   
  end if
  
  
  RETURN
  
END SUBROUTINE swap_faces
!
! NAME
!     SUBROUTINE iterative_improvement_points
!
! DESCRIPTION
!
!     
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE iterative_improvement_points(x1,y1,z1,x2,y2,z2,local_point_list,n_local_points)
  
USE local_mesh

IMPLICIT NONE
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  
  integer		:: n_local_points
  integer		:: local_point_list(1:n_local_points,1:3)

! local variables

  integer ix,iy,iz
  integer ixm,iym,izm
  integer ixp,iyp,izp
  integer ixo,iyo,izo
       
  real*8		:: dist_1,dist_2
  
  integer 		:: point
  integer 		:: n_moved
  
  integer,parameter	:: max_iterations=1000000
  integer		:: n_iterations
  
     
! START
	    
  n_iterations=0
  
!  write(*,*)'CALLED: iterative_improvement_points'

10  CONTINUE

    n_moved=0
    n_iterations=n_iterations+1
    
    if (n_iterations.GT.max_iterations) then
      write(*,*)'Maximum number of iterations exceeded in iterative_improvement_points'
      STOP
    end if
    
! loop over the line points looking for 'corner' points and 
! flip the corner point if this improves the fit between the meshed line and the original line
      
    do point=2,n_local_points-1

! check to see whether this is a corner point    

      ix=local_point_list(point,1)
      iy=local_point_list(point,2)
      iz=local_point_list(point,3)
      
      ixm=local_point_list(point-1,1)
      iym=local_point_list(point-1,2)
      izm=local_point_list(point-1,3)
      
      ixp=local_point_list(point+1,1)
      iyp=local_point_list(point+1,2)
      izp=local_point_list(point+1,3)
      
      CALL get_opposite_point(ixm,iym,izm,ix,iy,iz,ixp,iyp,izp,ixo,iyo,izo) 
      
      CALL dist_local_grid_point_to_line(ix ,iy ,iz ,x1,y1,z1,x2,y2,z2,dist_1)
      CALL dist_local_grid_point_to_line(ixo,iyo,izo,x1,y1,z1,x2,y2,z2,dist_2)
	     
      if (dist_2.LT.dist_1) then
        n_moved=n_moved+1
        local_point_list(point,1)=ixo
        local_point_list(point,2)=iyo
        local_point_list(point,3)=izo
!	write(*,*)n_moved,'Swap points',ix,iy,iz,ixo,iyo,izo
      end if
    
    end do
    
!    write(*,*)'n_moved=',n_moved

    if ((n_moved.NE.0)) GOTO 10
  
  RETURN
  
  END SUBROUTINE iterative_improvement_points
!
! NAME
!     SUBROUTINE check_edges
!
! DESCRIPTION
!     check_edges: the local_edge should be either 0 or 2 for a properly constructed triangle with no holes in it.
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/5/2013 CJS 
!
SUBROUTINE check_edges(error_found)
  
USE local_mesh

IMPLICIT NONE

  logical error_found

! local variables

  integer ix,iy,iz,value
   
! START

  error_found=.FALSE.

  do ix=local_ixmin,local_ixmax
    do iy=local_iymin,local_iymax
      do iz=local_izmin,local_izmax
      
        value=local_grid(ix,iy,iz,local_xedge)
	if ((value.ne.0).AND.(value.ne.2) ) then
!	  write(*,*)'Error found in check edges: x-edge value=',value
!	  write(*,*)'ix=',ix,' iy=',iy,' iz=',iz
          error_found=.TRUE.
	end if
      
        value=local_grid(ix,iy,iz,local_yedge)
	if ((value.ne.0).AND.(value.ne.2) ) then
!	  write(*,*)'Error found in check edges: y-edge value=',value
!	  write(*,*)'ix=',ix,' iy=',iy,' iz=',iz
          error_found=.TRUE.
	end if
      
        value=local_grid(ix,iy,iz,local_zedge)
	if ((value.ne.0).AND.(value.ne.2) ) then
!	  write(*,*)'Error found in check edges: z-edge value=',value
!	  write(*,*)'ix=',ix,' iy=',iy,' iz=',iz
          error_found=.TRUE.
	end if
      
      end do
    end do
  end do

  RETURN
  
  END SUBROUTINE check_edges
!
! NAME
!     SUBROUTINE add_edge
!
! DESCRIPTION
!     add an edge to the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE add_edge(ix1,iy1,iz1,ix2,iy2,iz2)
  
USE local_mesh

IMPLICIT NONE

  integer ix1,iy1,iz1,ix2,iy2,iz2

! local variables
  integer ox,oy,oz
  integer ix,iy,iz
  integer edge
   
! START

! add the face
!  write(*,*)'Add edge',ix1,iy1,iz1,ix2,iy2,iz2

  if      (ix2-ix1.eq.+1) then
    ix=ix1
    iy=iy1
    iz=iz1
    edge=local_xedge
  else if (ix2-ix1.eq.-1) then
    ix=ix2
    iy=iy1
    iz=iz1
    edge=local_xedge
  else if (iy2-iy1.eq.+1) then
    ix=ix1
    iy=iy1
    iz=iz1
    edge=local_yedge
  else if (iy2-iy1.eq.-1) then
    ix=ix1
    iy=iy2
    iz=iz1
    edge=local_yedge
  else if (iz2-iz1.eq.+1) then
    ix=ix1
    iy=iy1
    iz=iz1
    edge=local_zedge
  else if (iz2-iz1.eq.-1) then
    ix=ix1
    iy=iy1
    iz=iz2
    edge=local_zedge
  else
    write(*,*)'Error in add edge'
    write(*,*)'Not a valid edge'
    write(*,*)'point 1:',ix1,iy1,iz1
    write(*,*)'point 2:',ix2,iy2,iz2
    STOP
  end if
  
  local_grid(ix,iy,iz,edge)=local_grid(ix,iy,iz,edge)+1
  
  RETURN
  
END SUBROUTINE add_edge
!
! NAME
!     SUBROUTINE subtract_edge
!
! DESCRIPTION
!     subtract an edge to the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE subtract_edge(ix1,iy1,iz1,ix2,iy2,iz2)
  
USE local_mesh

IMPLICIT NONE

  integer ix1,iy1,iz1,ix2,iy2,iz2

! local variables
  integer ox,oy,oz
  integer ix,iy,iz
  integer edge
   
! START

! subtract the face
!  write(*,*)'subtract edge',ix1,iy1,iz1,ix2,iy2,iz2

  if      (ix2-ix1.eq.+1) then
    ix=ix1
    iy=iy1
    iz=iz1
    edge=local_xedge
  else if (ix2-ix1.eq.-1) then
    ix=ix2
    iy=iy1
    iz=iz1
    edge=local_xedge
  else if (iy2-iy1.eq.+1) then
    ix=ix1
    iy=iy1
    iz=iz1
    edge=local_yedge
  else if (iy2-iy1.eq.-1) then
    ix=ix1
    iy=iy2
    iz=iz1
    edge=local_yedge
  else if (iz2-iz1.eq.+1) then
    ix=ix1
    iy=iy1
    iz=iz1
    edge=local_zedge
  else if (iz2-iz1.eq.-1) then
    ix=ix1
    iy=iy1
    iz=iz2
    edge=local_zedge
  else
    write(*,*)'Error in subtract edge'
    write(*,*)'Not a valid edge'
    write(*,*)'point 1:',ix1,iy1,iz1
    write(*,*)'point 2:',ix2,iy2,iz2
    STOP
  end if
  
  local_grid(ix,iy,iz,edge)=local_grid(ix,iy,iz,edge)-1
 
  RETURN
   
END SUBROUTINE subtract_edge
!
! NAME
!     SUBROUTINE get_opposite_point
!
! DESCRIPTION
!     subtract an edge to the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE get_opposite_point(ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3,ox,oy,oz)
  
USE local_mesh

IMPLICIT NONE

  integer ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3,ox,oy,oz

! local variables
  integer offset1_x,offset1_y,offset1_z
  integer offset2_x,offset2_y,offset2_z
   
! START

  offset1_x=ix2-ix1
  offset1_y=iy2-iy1
  offset1_z=iz2-iz1

  offset2_x=ix3-ix2
  offset2_y=iy3-iy2
  offset2_z=iz3-iz2
  
  ox=ix1+offset2_x
  oy=iy1+offset2_y
  oz=iz1+offset2_z
    
  RETURN
  
END SUBROUTINE get_opposite_point
!
! NAME
!     SUBROUTINE swap_points
!
! DESCRIPTION
!     subtract an edge to the local_grid
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE swap_points(ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3,ox,oy,oz)
  
USE local_mesh

IMPLICIT NONE

  integer ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3,ox,oy,oz

! local variables
   
! START
  
  local_grid(ix2,iy2,iz2,local_corner)=0
!  CALL subtract_edge(ix1,iy1,iz1,ix2,iy2,iz2)
!  CALL subtract_edge(ix2,iy2,iz2,ix3,iy3,iz3)
 
  local_grid(ox ,oy ,oz ,local_corner)=1
!  CALL add_edge(ix1,iy1,iz1,ox,oy,oz)
!  CALL add_edge(ox,oy,oz,ix3,iy3,iz3)
    
  RETURN
  
END SUBROUTINE swap_points
!
! NAME
!     SUBROUTINE order_points_2
!
! DESCRIPTION
!     order_points_2:
!
!     put two points in order defined by their x, y and z coordinates
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE order_points_2(x1_in,y1_in,z1_in,x2_in,y2_in,z2_in,x1,y1,z1,x2,y2,z2)
  
USE local_mesh

IMPLICIT NONE
  
  real*8		:: x1_in,y1_in,z1_in
  real*8		:: x2_in,y2_in,z2_in

  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2

! local variables

  integer		:: point_order
  
! START

! work out the ordering of the points i.e. which way to traverse a line defined by two points
! Ordering the points should make the mesh generation more consistent

  point_order=1
  
! look to traverse in +x direction first  
  if (x1_in.lt.x2_in) then
    point_order=1
  else if (x1_in.gt.x2_in) then
    point_order=2
  else
  
! look to traverse in +y direction  
    if (y1_in.lt.y2_in) then
      point_order=1
    else if (y1_in.gt.y2_in) then
      point_order=2
    else
    
! look to traverse in +z direction  
      if (z1_in.lt.z2_in) then
        point_order=1
      else if (z1_in.gt.z2_in) then
        point_order=2
      else
! the end points must be the same
        point_order=1      
      end if ! look to traverse in +z direction 
      
    end if ! look to traverse in +y direction  

  end if ! look to traverse in +x direction 

  if (point_order.eq.1) then
    x1=x1_in
    y1=y1_in
    z1=z1_in
    x2=x2_in
    y2=y2_in
    z2=z2_in
  else
    x1=x2_in
    y1=y2_in
    z1=z2_in
    x2=x1_in
    y2=y1_in
    z2=z1_in
  end if
   
  RETURN
    
  END SUBROUTINE order_points_2
!
! NAME
!     SUBROUTINE order_points_2
!
! DESCRIPTION
!     order_points_2:
!
!     put two points in order defined by their x, y and z coordinates
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 3/6/2013 CJS 
!
SUBROUTINE order_points_3(x1_in,y1_in,z1_in,x2_in,y2_in,z2_in,x3_in,y3_in,z3_in,	&
                          x1,y1,z1,x2,y2,z2,x3,y3,z3)
  
USE local_mesh

IMPLICIT NONE
  
  real*8		:: x1_in,y1_in,z1_in
  real*8		:: x2_in,y2_in,z2_in
  real*8		:: x3_in,y3_in,z3_in

  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3

! local variables
  
! START

! work out the ordering of the points i.e. which way to traverse a triangle defined by three points
! Ordering the points should make the mesh generation more consistent

! copy points
  x1=x1_in
  y1=y1_in
  z1=z1_in
  
  x2=x2_in
  y2=y2_in
  z2=z2_in
  
  x3=x3_in
  y3=y3_in
  z3=z3_in
  
  CALL order_points_2(x2,y2,z2,x3,y3,z3,x2,y2,z2,x3,y3,z3)
  CALL order_points_2(x1,y1,z1,x2,y2,z2,x1,y1,z1,x2,y2,z2)
  CALL order_points_2(x2,y2,z2,x3,y3,z3,x2,y2,z2,x3,y3,z3)
   
  RETURN
    
  END SUBROUTINE order_points_3
