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
!SUBROUTINE build_surface_sphere
!SUBROUTINE scale_point
!SUBROUTINE get_point_from_list
!SUBROUTINE build_surface_sphere_OLD
!SUBROUTINE build_surface_sphere_OLD_OLD
!

! NAME
!     SUBROUTINE build_surface_sphere
!
! DESCRIPTION
!     build_surface_sphere:
!
!     create a triangulated sphere surface
!     
! COMMENTS
!   
! Edges and points associated with triangles should be defined in a clockwise direction
! around the triangle at each stage. The sign of an edge in the triangle_to_edge data structure
! indicates the direction in which an edge should be traversed.
!
! The process works using triangle subdivision with the starting point of an icosahedron
! note that the icosohedron altitude angle is not exact - this needs to be worked out properly.
! See http://www.vb-helper.com/tutorial_platonic_solids.html for example...
!
! HISTORY
!
!     started 29/08/2012 CJS
!     surface roughness 13/8/2015 CJS
!     use triangle subdivision method 21/9/2015 CJS
!     Fix error which in which transformations (translations and rotations were in error) CJS 12/11/2019
!
!
SUBROUTINE build_surface_sphere(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE geometry_operators
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: theta,phi
real*8	:: dtheta,dphi
real*8	:: radius0,radius

real*8	:: edge_shrink_factor,r_n_stages

integer	:: stage,n_stages

integer		:: initial_n_points,new_n_points,tot_n_points
integer		:: point_count,last_point_count,point

integer		:: initial_n_triangles,new_n_triangles,tot_n_triangles
integer		:: triangle_count,last_triangle_count,triangle
	
integer		:: initial_n_edges,new_n_edges,tot_n_edges
integer		:: edge_count,last_edge_count,edge,new_edge
integer		:: e1,e2,e3
integer		:: elist(9)
integer		:: sign

type(xyz)	:: point1,point2,point3,point4,point5,point6
integer		:: p1,p2,p3,p4,p5,p6
integer		:: p1check,p2check,p3check,p4check,p5check,p6check

type(xyz),allocatable	:: point_list(:)

integer,allocatable	:: triangle_to_point_list(:,:)
integer,allocatable	:: new_triangle_to_point_list(:,:)

integer,allocatable	:: triangle_to_edge_list(:,:)
integer,allocatable	:: new_triangle_to_edge_list(:,:)

integer,allocatable	:: edge_to_point_list(:,:)
integer,allocatable	:: new_edge_to_point_list(:,:)

real*8		:: r_random
integer		:: tvalue
integer		:: pvalue
integer		:: i,ecount
integer		:: i1,i2

integer		:: check_OK
logical		:: check_failed

! START

    
! CREATE A TRIANGULATED SPHERE MESH ON THE SCALE OF DL.    
      
      radius0=problem_surfaces(surface_number)%surface_parameters(1)

! The length of triangle edges to be of the order of dl/2, half the mesh edge length 

! calculate the number of iterative triangle sub-divisions required to achieve the final edge length
      edge_shrink_factor=radius0/(dl/2d0)
      r_n_stages=log(edge_shrink_factor)/log(2d0)
      n_stages=NINT(r_n_stages)
      if (n_stages.LE.1) n_stages=2
      
      write(*,*)'radius=',radius0,' dl/2=',dl/2d0
      write(*,*)'Number of stages=',n_stages
      
! calculate the number of points, triangles and edges and allocate memory for the point data

! we start with an icosahedron i.e. 6 points, 8 triangles and 12 edges
! then at each subdivision stage:
! 1: each edge gets split into 2 to introduce a new point
! 2: each triangle get subdivided into 4 i.e. 3 new ones are added
! 3: each edge gets split into two and three new edges are added for each oof the original triangles

      initial_n_points=12
      initial_n_triangles=20   
      initial_n_edges=30
     
      do stage=1,n_stages
      
        new_n_points    =initial_n_points+initial_n_edges
        new_n_triangles =initial_n_triangles+3*initial_n_triangles
        new_n_edges     =initial_n_edges+initial_n_edges+3*initial_n_triangles

        initial_n_points   =new_n_points
        initial_n_triangles=new_n_triangles
        initial_n_edges    =new_n_edges
		
      end do ! next stage
      
      tot_n_points    =initial_n_points
      tot_n_triangles =initial_n_triangles
      tot_n_edges     =initial_n_edges
      
      write(*,*)'tot_n_points',tot_n_points
      write(*,*)'tot_n_triangles',tot_n_triangles
      write(*,*)'tot_n_edges',tot_n_edges
	
      problem_surfaces(surface_number)%number_of_triangles=tot_n_triangles 
      
      ALLOCATE( point_list(1:tot_n_points) )
      ALLOCATE( triangle_to_point_list(1:tot_n_triangles,1:3) )
      ALLOCATE( new_triangle_to_point_list(1:tot_n_triangles,1:3) )
      ALLOCATE( triangle_to_edge_list(1:tot_n_triangles,1:3) )
      ALLOCATE( new_triangle_to_edge_list(1:tot_n_triangles,1:3) )
      ALLOCATE( edge_to_point_list(1:tot_n_edges,1:4) )
      ALLOCATE( new_edge_to_point_list(1:tot_n_edges,1:4) )
      
      ALLOCATE( problem_surfaces(surface_number)%triangle_list(1:tot_n_triangles) )
      
! calculate the initial points and triangles
! north pole point      
      point_count=1
      theta=0d0
      phi=0d0
      
      radius=radius0
      if (problem_surfaces(surface_number)%roughness_flag) then
        CALL random_number(r_random)
	r_random=(r_random-0.5d0)*2d0
        radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
      end if
      CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
!      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      point_list(point_count)=point1
      
! south pole point      
      point_count=point_count+1
      theta=pi
      phi=0d0
      
      radius=radius0
      if (problem_surfaces(surface_number)%roughness_flag) then
        CALL random_number(r_random)
	r_random=(r_random-0.5d0)*2d0
        radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
      end if
      CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
!      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      point_list(point_count)=point1
  
! top layer of 5 points     
      theta= pi/3d0   !****** THIS ANGLE IS ONLY APPROXIMATE - TO BE DEFINED PROPERLY******
      do i=1,5
      
        point_count=point_count+1
        phi=(i-1)*pi*2d0/5d0
     
        radius=radius0
        if (problem_surfaces(surface_number)%roughness_flag) then
     	  CALL random_number(r_random)
	  r_random=(r_random-0.5d0)*2d0
     	  radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
        end if
        CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
!        CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
        point_list(point_count)=point1
      
      end do ! next equatorial point
  
! bottom layer of 5 points     
      theta= 2d0*pi/3d0   !****** THIS ANGLE IS ONLY APPROXIMATE - TO BE DEFINED PROPERLY******
      do i=1,5
      
        point_count=point_count+1
        phi=(i-1)*pi*2d0/5d0+pi/5d0
     
        radius=radius0
        if (problem_surfaces(surface_number)%roughness_flag) then
     	  CALL random_number(r_random)
	  r_random=(r_random-0.5d0)*2d0
     	  radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
        end if
        CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
!        CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
        point_list(point_count)=point1
      
      end do ! next equatorial point

! set the triangle to point list and the triangle to edge list
      triangle_to_point_list(1:tot_n_triangles,1:3)=0
      triangle_to_edge_list(1:tot_n_triangles,1:3)=0
      
! set triangles to make the normal point outwards, note a negative edge number indicates that it is traversed backwards	  
      triangle_count=1
      triangle_to_point_list(triangle_count,1)=1
      triangle_to_point_list(triangle_count,2)=3
      triangle_to_point_list(triangle_count,3)=4
      
      triangle_to_edge_list(triangle_count,1)= 1
      triangle_to_edge_list(triangle_count,2)= 6
      triangle_to_edge_list(triangle_count,3)=-2
      
! next triangle  (2)
      triangle_count=2
      triangle_to_point_list(triangle_count,1)=1
      triangle_to_point_list(triangle_count,2)=4
      triangle_to_point_list(triangle_count,3)=5
      
      triangle_to_edge_list(triangle_count,1)= 2
      triangle_to_edge_list(triangle_count,2)= 7
      triangle_to_edge_list(triangle_count,3)=-3
      
! next triangle      (3)
      triangle_count=3
      triangle_to_point_list(triangle_count,1)=1
      triangle_to_point_list(triangle_count,2)=5
      triangle_to_point_list(triangle_count,3)=6
      
      triangle_to_edge_list(triangle_count,1)= 3
      triangle_to_edge_list(triangle_count,2)= 8
      triangle_to_edge_list(triangle_count,3)=-4
      
! next triangle      (4)
      triangle_count=4
      triangle_to_point_list(triangle_count,1)=1
      triangle_to_point_list(triangle_count,2)=6
      triangle_to_point_list(triangle_count,3)=7
      
      triangle_to_edge_list(triangle_count,1)= 4
      triangle_to_edge_list(triangle_count,2)= 9
      triangle_to_edge_list(triangle_count,3)=-5
      
! next triangle      (5)
      triangle_count=5
      triangle_to_point_list(triangle_count,1)=1
      triangle_to_point_list(triangle_count,2)=7
      triangle_to_point_list(triangle_count,3)=3
      
      triangle_to_edge_list(triangle_count,1)= 5
      triangle_to_edge_list(triangle_count,2)= 10
      triangle_to_edge_list(triangle_count,3)=-1
      
! next triangle      (6)
      triangle_count=6
      triangle_to_point_list(triangle_count,1)=3
      triangle_to_point_list(triangle_count,2)=8
      triangle_to_point_list(triangle_count,3)=4
      
      triangle_to_edge_list(triangle_count,1)= 11
      triangle_to_edge_list(triangle_count,2)=-12
      triangle_to_edge_list(triangle_count,3)=-6
      
! next triangle      (7)
      triangle_count=7
      triangle_to_point_list(triangle_count,1)=4
      triangle_to_point_list(triangle_count,2)=8
      triangle_to_point_list(triangle_count,3)=9
      
      triangle_to_edge_list(triangle_count,1)= 12
      triangle_to_edge_list(triangle_count,2)= 21
      triangle_to_edge_list(triangle_count,3)=-13
      
! next triangle      (8)
      triangle_count=8
      triangle_to_point_list(triangle_count,1)=5
      triangle_to_point_list(triangle_count,2)=4
      triangle_to_point_list(triangle_count,3)=9
      
      triangle_to_edge_list(triangle_count,1)=-7
      triangle_to_edge_list(triangle_count,2)= 13
      triangle_to_edge_list(triangle_count,3)=-14
      
! next triangle      (9)
      triangle_count=9
      triangle_to_point_list(triangle_count,1)=5
      triangle_to_point_list(triangle_count,2)=9
      triangle_to_point_list(triangle_count,3)=10
      
      triangle_to_edge_list(triangle_count,1)= 14
      triangle_to_edge_list(triangle_count,2)= 22
      triangle_to_edge_list(triangle_count,3)=-15
      
! next triangle      (10)
      triangle_count=10
      triangle_to_point_list(triangle_count,1)=6
      triangle_to_point_list(triangle_count,2)=5
      triangle_to_point_list(triangle_count,3)=10
      
      triangle_to_edge_list(triangle_count,1)=-8
      triangle_to_edge_list(triangle_count,2)= 15
      triangle_to_edge_list(triangle_count,3)=-16
      
! next triangle      (11)
      triangle_count=11
      triangle_to_point_list(triangle_count,1)=6
      triangle_to_point_list(triangle_count,2)=10
      triangle_to_point_list(triangle_count,3)=11
      
      triangle_to_edge_list(triangle_count,1)= 16
      triangle_to_edge_list(triangle_count,2)= 23
      triangle_to_edge_list(triangle_count,3)=-17
      
! next triangle      (12)
      triangle_count=12
      triangle_to_point_list(triangle_count,1)=7
      triangle_to_point_list(triangle_count,2)=6
      triangle_to_point_list(triangle_count,3)=11
      
      triangle_to_edge_list(triangle_count,1)=-9
      triangle_to_edge_list(triangle_count,2)= 17
      triangle_to_edge_list(triangle_count,3)=-18
      
! next triangle      (13)
      triangle_count=13
      triangle_to_point_list(triangle_count,1)=7
      triangle_to_point_list(triangle_count,2)=11
      triangle_to_point_list(triangle_count,3)=12
      
      triangle_to_edge_list(triangle_count,1)= 18
      triangle_to_edge_list(triangle_count,2)= 24
      triangle_to_edge_list(triangle_count,3)=-19
      
! next triangle      (14)
      triangle_count=14
      triangle_to_point_list(triangle_count,1)=3
      triangle_to_point_list(triangle_count,2)=7
      triangle_to_point_list(triangle_count,3)=12
      
      triangle_to_edge_list(triangle_count,1)=-10
      triangle_to_edge_list(triangle_count,2)= 19
      triangle_to_edge_list(triangle_count,3)=-20
      
! next triangle      (15)
      triangle_count=15
      triangle_to_point_list(triangle_count,1)=3
      triangle_to_point_list(triangle_count,2)=12
      triangle_to_point_list(triangle_count,3)=8
      
      triangle_to_edge_list(triangle_count,1)= 20
      triangle_to_edge_list(triangle_count,2)= 25
      triangle_to_edge_list(triangle_count,3)=-11
      
! next triangle      (16)
      triangle_count=16
      triangle_to_point_list(triangle_count,1)=9
      triangle_to_point_list(triangle_count,2)=8
      triangle_to_point_list(triangle_count,3)=2
      
      triangle_to_edge_list(triangle_count,1)=-21
      triangle_to_edge_list(triangle_count,2)= 26
      triangle_to_edge_list(triangle_count,3)=-27
      
! next triangle      (17)
      triangle_count=17
      triangle_to_point_list(triangle_count,1)=10
      triangle_to_point_list(triangle_count,2)=9
      triangle_to_point_list(triangle_count,3)=2
      
      triangle_to_edge_list(triangle_count,1)=-22
      triangle_to_edge_list(triangle_count,2)= 27
      triangle_to_edge_list(triangle_count,3)=-28
      
! next triangle      (18)
      triangle_count=18
      triangle_to_point_list(triangle_count,1)=11
      triangle_to_point_list(triangle_count,2)=10
      triangle_to_point_list(triangle_count,3)=2
      
      triangle_to_edge_list(triangle_count,1)=-23
      triangle_to_edge_list(triangle_count,2)= 28
      triangle_to_edge_list(triangle_count,3)=-29
      
! next triangle      (19)
      triangle_count=19
      triangle_to_point_list(triangle_count,1)=12
      triangle_to_point_list(triangle_count,2)=11
      triangle_to_point_list(triangle_count,3)=2
      
      triangle_to_edge_list(triangle_count,1)=-24
      triangle_to_edge_list(triangle_count,2)= 29
      triangle_to_edge_list(triangle_count,3)=-30
      
! next triangle      (20)
      triangle_count=20
      triangle_to_point_list(triangle_count,1)=8
      triangle_to_point_list(triangle_count,2)=12
      triangle_to_point_list(triangle_count,3)=2
      
      triangle_to_edge_list(triangle_count,1)=-25
      triangle_to_edge_list(triangle_count,2)= 30
      triangle_to_edge_list(triangle_count,3)=-26
      
! set the edge to point list
      edge_to_point_list(1:tot_n_edges,1:4)=0
      
      edge_count=1
      edge_to_point_list(edge_count,1)=1
      edge_to_point_list(edge_count,2)=3
      
      edge_count=2
      edge_to_point_list(edge_count,1)=1
      edge_to_point_list(edge_count,2)=4
      
      edge_count=3
      edge_to_point_list(edge_count,1)=1
      edge_to_point_list(edge_count,2)=5
      
      edge_count=4
      edge_to_point_list(edge_count,1)=1
      edge_to_point_list(edge_count,2)=6
      
      edge_count=5
      edge_to_point_list(edge_count,1)=1
      edge_to_point_list(edge_count,2)=7
      
      edge_count=6
      edge_to_point_list(edge_count,1)=3
      edge_to_point_list(edge_count,2)=4
      
      edge_count=7
      edge_to_point_list(edge_count,1)=4
      edge_to_point_list(edge_count,2)=5
      
      edge_count=8
      edge_to_point_list(edge_count,1)=5
      edge_to_point_list(edge_count,2)=6
      
      edge_count=9
      edge_to_point_list(edge_count,1)=6
      edge_to_point_list(edge_count,2)=7
      
      edge_count=10
      edge_to_point_list(edge_count,1)=7
      edge_to_point_list(edge_count,2)=3
      
      edge_count=11
      edge_to_point_list(edge_count,1)=3
      edge_to_point_list(edge_count,2)=8
       
      edge_count=12
      edge_to_point_list(edge_count,1)=4
      edge_to_point_list(edge_count,2)=8
      
      edge_count=13
      edge_to_point_list(edge_count,1)=4
      edge_to_point_list(edge_count,2)=9
       
      edge_count=14
      edge_to_point_list(edge_count,1)=5
      edge_to_point_list(edge_count,2)=9
      
      edge_count=15
      edge_to_point_list(edge_count,1)=5
      edge_to_point_list(edge_count,2)=10
       
      edge_count=16
      edge_to_point_list(edge_count,1)=6
      edge_to_point_list(edge_count,2)=10
      
      edge_count=17
      edge_to_point_list(edge_count,1)=6
      edge_to_point_list(edge_count,2)=11
       
      edge_count=18
      edge_to_point_list(edge_count,1)=7
      edge_to_point_list(edge_count,2)=11
      
      edge_count=19
      edge_to_point_list(edge_count,1)=7
      edge_to_point_list(edge_count,2)=12
       
      edge_count=20
      edge_to_point_list(edge_count,1)=3
      edge_to_point_list(edge_count,2)=12
       
      edge_count=21
      edge_to_point_list(edge_count,1)=8
      edge_to_point_list(edge_count,2)=9
       
      edge_count=22
      edge_to_point_list(edge_count,1)=9
      edge_to_point_list(edge_count,2)=10
       
      edge_count=23
      edge_to_point_list(edge_count,1)=10
      edge_to_point_list(edge_count,2)=11
       
      edge_count=24
      edge_to_point_list(edge_count,1)=11
      edge_to_point_list(edge_count,2)=12
       
      edge_count=25
      edge_to_point_list(edge_count,1)=12
      edge_to_point_list(edge_count,2)=8
       
      edge_count=26
      edge_to_point_list(edge_count,1)=8
      edge_to_point_list(edge_count,2)=2
       
      edge_count=27
      edge_to_point_list(edge_count,1)=9
      edge_to_point_list(edge_count,2)=2
       
      edge_count=28
      edge_to_point_list(edge_count,1)=10
      edge_to_point_list(edge_count,2)=2
       
      edge_count=29
      edge_to_point_list(edge_count,1)=11
      edge_to_point_list(edge_count,2)=2
       
      edge_count=30
      edge_to_point_list(edge_count,1)=12
      edge_to_point_list(edge_count,2)=2
     

! LOOP OVER ITERATIVE TRIANGLE SUB-DIVISIONS TO CREATE NEW TRIANGLES
      
      do stage=1,n_stages
      
	check_failed=.FALSE.   ! reset edge check for this stage of subdivision
	
        new_triangle_to_point_list(1:tot_n_triangles,1:3)=0
        new_triangle_to_edge_list(1:tot_n_triangles,1:3)=0
        new_edge_to_point_list(1:tot_n_edges,1:4)=0    

! copy and reset counters for new triangles, edges and points	

        last_triangle_count=triangle_count
        last_edge_count=edge_count
        last_point_count=point_count
        triangle_count=0
        edge_count=0
! note don't reset point_count as we just add to the original array
      
! PROCESS 1: EACH EDGE GETS SPLIT INTO 2 TO INTRODUCE A NEW POINT
      
        do edge=1,last_edge_count
	
! get the two points on the edge

          p1=edge_to_point_list(edge,1)
          p2=edge_to_point_list(edge,2)
	  point1=point_list(p1)
	  point2=point_list(p2)
	  
! calculate the new point at the edge centre, projected out onto the sphere
          point_count=point_count+1
	  p3=point_count
	  point3=point1+point2
! scale to radius
          if (problem_surfaces(surface_number)%roughness_flag) then
            CALL random_number(r_random)
	    r_random=(r_random-0.5d0)*2d0
            radius=radius0+r_random*problem_surfaces(surface_number)%roughness_p1
          end if
	  
          CALL scale_point(radius,point3)
	  
!          CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
	  point_list(p3)=point3
	  edge_to_point_list(edge,3)=p3
	  
	end do ! next edge
	
! PROCESS 2: EACH TRIANGLE GET SUBDIVIDED INTO 4 I.E. 3 NEW ONES ARE ADDED 
	
        do triangle=1,last_triangle_count

! 2A: get the original 3 points on this triangle
          p1=triangle_to_point_list(triangle,1)
          p2=triangle_to_point_list(triangle,2)
          p3=triangle_to_point_list(triangle,3)
	   
! 2B: get the 3 new points on this triangle as the new points splitting the triangle edges
! note: use abs(e) as the edge number can be negative if it is traversed in the opposite direction to the point definition

          e1=abs(triangle_to_edge_list(triangle,1))
	  p4=edge_to_point_list(e1,3)
          e2=abs(triangle_to_edge_list(triangle,2))
	  p5=edge_to_point_list(e2,3)
          e3=abs(triangle_to_edge_list(triangle,3))
	  p6=edge_to_point_list(e3,3)
	  
! Check that the trinagle edges and points are oriented correctly...

          check_OK=0

! edge 1 check
	  if (triangle_to_edge_list(triangle,1).LT.0) then
	    p1check=edge_to_point_list(e1,2)
	    p2check=edge_to_point_list(e1,1)
	  else
	    p1check=edge_to_point_list(e1,1)
	    p2check=edge_to_point_list(e1,2)
	  end if
	  
	  if ( (p1.NE.p1check).OR.(p2.NE.p2check) ) check_OK=1

! edge 2 check
	  if (triangle_to_edge_list(triangle,2).LT.0) then
	    p2check=edge_to_point_list(e2,2)
	    p3check=edge_to_point_list(e2,1)
	  else
	    p2check=edge_to_point_list(e2,1)
	    p3check=edge_to_point_list(e2,2)
	  end if
	  
	  if ( (p2.NE.p2check).OR.(p3.NE.p3check) ) check_OK=2

! edge 3 check
	  if (triangle_to_edge_list(triangle,3).LT.0) then
	    p3check=edge_to_point_list(e3,2)
	    p1check=edge_to_point_list(e3,1)
	  else
	    p3check=edge_to_point_list(e3,1)
	    p1check=edge_to_point_list(e3,2)
	  end if
	  
	  if ( (p3.NE.p3check).OR.(p1.NE.p1check) ) check_OK=3
	  
	  if (check_OK.NE.0) then
	    write(*,*)'Failed edge and point consistency check on edge:',check_OK
	    write(*,*)'Stage:',stage
	    write(*,*)'Triangle number',triangle
	    
	    write(*,*)'Point numbers    ',p1,p2,p3
	    write(*,*)'Edge 1',triangle_to_edge_list(triangle,1)
	    write(*,*)'End point numbers',edge_to_point_list(e1,1),edge_to_point_list(e1,2)
	    write(*,*)'Edge 2',triangle_to_edge_list(triangle,2)
	    write(*,*)'End point numbers',edge_to_point_list(e2,1),edge_to_point_list(e2,2)
	    write(*,*)'Edge 3',triangle_to_edge_list(triangle,3)
	    write(*,*)'End point numbers',edge_to_point_list(e3,1),edge_to_point_list(e3,2)
	    check_failed=.TRUE.
	  end if  

! 2C: loop over the original edges in this triangle, splitting them if required to create the new edge list

          ecount=0
          do i=1,3
	  
	    edge=abs(triangle_to_edge_list(triangle,i))
	    if (triangle_to_edge_list(triangle,i).LT.0) then
	      sign=-1
	    else
	      sign=1
	    end if
	    
! check whether this edge has been split before, if not split it now.
! the new edge is placed in element 4 of the edge_to_point_list array

            if (edge_to_point_list(edge,4).EQ.0) then
	    
! this edge is not split so create two edges in the new_edge_to_point_list structure
! first edge
              edge_count=edge_count+1
	      e1=edge_count
	      new_edge_to_point_list(edge_count,1)=edge_to_point_list(edge,1) ! first point
	      new_edge_to_point_list(edge_count,2)=edge_to_point_list(edge,3) ! middle point
	      new_edge_to_point_list(edge_count,3)=0
	      new_edge_to_point_list(edge_count,4)=0

! second edge
              edge_count=edge_count+1
	      e2=edge_count
	      new_edge_to_point_list(edge_count,1)=edge_to_point_list(edge,3) ! middle point
	      new_edge_to_point_list(edge_count,2)=edge_to_point_list(edge,2) ! second point
	      new_edge_to_point_list(edge_count,3)=0
	      new_edge_to_point_list(edge_count,4)=0

! set the first new edge number in the original edge list	
       
	      edge_to_point_list(edge,4)=e1

	    end if
	    
	    if (sign.eq.1) then
	      e1=edge_to_point_list(edge,4)
	      e2=e1+1
	    else
	      e2=edge_to_point_list(edge,4)
	      e1=e2+1	
	    end if

! add the edges to the local list of edges, taking into account the direction of the edge
	    ecount=ecount+1
	    elist(ecount)=sign*e1
	    
	    ecount=ecount+1
	    elist(ecount)=sign*e2
	    
	  end do ! next original edge of the triangle

! create three new internal edges in this triangle
	  
          edge_count=edge_count+1
! set points on the new edge	      
	  new_edge_to_point_list(edge_count,1)=p4
	  new_edge_to_point_list(edge_count,2)=p5
	  
! add the edges to the local list of edges, taking into account the direction of the edge
	  ecount=ecount+1
	  elist(ecount)=edge_count
	  
          edge_count=edge_count+1
! set points on the new edge	      
	  new_edge_to_point_list(edge_count,1)=p5
	  new_edge_to_point_list(edge_count,2)=p6
	  
! add the edges to the local list of edges, taking into account the direction of the edge
	  ecount=ecount+1
	  elist(ecount)=edge_count
	  
          edge_count=edge_count+1
! set points on the new edge	      
	  new_edge_to_point_list(edge_count,1)=p6
	  new_edge_to_point_list(edge_count,2)=p4
	  
! add the edges to the local list of edges, taking into account the direction of the edge
	  ecount=ecount+1
	  elist(ecount)=edge_count
	  
! at this stage we have the point list p1-p6
! and the edge list elist(1:9)
! we can now create the new triangles and the new triangle to point list and triangle to edge list
		  
! generate the new triangles     

          triangle_count=triangle_count+1
          new_triangle_to_point_list(triangle_count,1)=p1
          new_triangle_to_point_list(triangle_count,2)=p4
          new_triangle_to_point_list(triangle_count,3)=p6
          new_triangle_to_edge_list(triangle_count,1)= elist(1)
          new_triangle_to_edge_list(triangle_count,2)=-elist(9)
          new_triangle_to_edge_list(triangle_count,3)= elist(6)

          triangle_count=triangle_count+1
	  new_triangle_to_point_list(triangle_count,1)=p6
          new_triangle_to_point_list(triangle_count,2)=p5
          new_triangle_to_point_list(triangle_count,3)=p3
          new_triangle_to_edge_list(triangle_count,1)=-elist(8)
          new_triangle_to_edge_list(triangle_count,2)= elist(4)
          new_triangle_to_edge_list(triangle_count,3)= elist(5)

          triangle_count=triangle_count+1
	  new_triangle_to_point_list(triangle_count,1)=p6
          new_triangle_to_point_list(triangle_count,2)=p4
          new_triangle_to_point_list(triangle_count,3)=p5
          new_triangle_to_edge_list(triangle_count,1)= elist(9)
          new_triangle_to_edge_list(triangle_count,2)= elist(7)
          new_triangle_to_edge_list(triangle_count,3)= elist(8)

          triangle_count=triangle_count+1
	  new_triangle_to_point_list(triangle_count,1)=p4
          new_triangle_to_point_list(triangle_count,2)=p2
          new_triangle_to_point_list(triangle_count,3)=p5
          new_triangle_to_edge_list(triangle_count,1)= elist(2)
          new_triangle_to_edge_list(triangle_count,2)= elist(3)
          new_triangle_to_edge_list(triangle_count,3)=-elist(7)
	  
        end do ! next triangle
	
! set up for the next subdivision stage	
        last_point_count=point_count
        last_triangle_count=triangle_count
        last_edge_count=edge_count
	
	triangle_to_point_list(1:tot_n_triangles,1:3)=new_triangle_to_point_list(1:tot_n_triangles,1:3)
	
	triangle_to_edge_list(1:tot_n_triangles,1:3)=new_triangle_to_edge_list(1:tot_n_triangles,1:3)
	
        edge_to_point_list(1:tot_n_edges,1:2)=new_edge_to_point_list(1:tot_n_edges,1:2)
        edge_to_point_list(1:tot_n_edges,3:4)=0
  
        if (check_failed) then
          write(*,*)'**** EDGE CHECK FAILED ****'
	  STOP
        end if
	  
      end do ! next stage of sub-division
      
! Apply the transformation to all the points in the point list

      do i=1,tot_n_points
      
        point1=point_list(i)
        CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
        point_list(i)=point1
       
      end do

! Copy the triangulated surface into the problem_surfaces structure      
      do triangle=1,triangle_count
        do i=1,3
	  problem_surfaces(surface_number)%triangle_list(triangle)%vertex(i)=point_list(triangle_to_point_list(triangle,i))
	end do
      end do

! Deallocate local memory
      DEALLOCATE( point_list )
      DEALLOCATE( triangle_to_point_list )
      DEALLOCATE( new_triangle_to_point_list )
      DEALLOCATE( triangle_to_edge_list )
      DEALLOCATE( new_triangle_to_edge_list )
      DEALLOCATE( edge_to_point_list )
      DEALLOCATE( new_edge_to_point_list )

  RETURN

END SUBROUTINE build_surface_sphere

!
! NAME
!     SUBROUTINE scale_point
!
! DESCRIPTION
!    a point is given which defines a vector from the origin. Scale the 
!    point position such that the vector is a given length
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 21/9/2015 CJS
!
!
SUBROUTINE scale_point(length,point)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

real*8		:: length
type(xyz)	:: point


! local variables

real*8		:: original_length

! START
		
  original_length=sqrt(point%x**2+point%y**2+point%z**2)
  point%x=point%x*length/original_length
  point%y=point%y*length/original_length
  point%z=point%z*length/original_length
  
END SUBROUTINE scale_point
!
! NAME
!     SUBROUTINE get_point_from_list
!
! DESCRIPTION
!    get_point_from_list
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/8/2015 CJS
!
!
SUBROUTINE get_point_from_list(theta,phi,n_theta,n_phi,point,point_list,tot_n_points)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer :: theta,phi,n_theta,n_phi

type(xyz)	:: point

integer		:: tot_n_points
type(xyz)	:: point_list(1:tot_n_points)


! local variables

integer :: lphi,lpoint

! START
	
  if (theta.eq.1) then
  
! point is at the north pole
    point=point_list(1)
    
  else if (theta.eq.n_theta+1) then
  
! point is at the south pole
    point=point_list(tot_n_points)
    
  else
  
    if (phi.EQ.n_phi+1) then	    
      lphi=1	      
    else 	    
      lphi=phi	      
    end if
    
    lpoint=2+(theta-2)*n_phi+(lphi-1)
    point=point_list(lpoint)
    
  end if
  
END SUBROUTINE get_point_from_list
!
! NAME
!     SUBROUTINE build_surface_sphere
!
! DESCRIPTION
!     build_surface_sphere_OLD:
!
!     create a triangulated sphere surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!     surface roughness 13/8/2015 CJS
!
!
SUBROUTINE build_surface_sphere_OLD(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: theta,phi
real*8	:: dtheta,dphi
real*8	:: radius0,radius
real*8	:: circumference
integer	:: n_circumference
integer :: n_theta,n_phi
integer :: theta_loop,phi_loop

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3,point4

integer		:: tot_n_points
integer		:: point_count
type(xyz),allocatable	:: point_list(:)
real*8		:: r_random
integer		:: tvalue
integer		:: pvalue

! START

    
! CREATE A TRIANGULATED SPHERE MESH ON THE SCALE OF DL.    
      
      radius0=problem_surfaces(surface_number)%surface_parameters(1)

! set the length of triangle edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      circumference=2d0*pi*radius0
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_theta=n_circumference/2
      n_phi=n_circumference
      
      dtheta=pi/n_theta
      dphi=2d0*pi/n_phi
      
! calculate the number of points and allocate memory for the point data
      tot_n_points=1  				! north pole point 
      tot_n_points=tot_n_points+(n_theta-1)*n_phi  	! surface points
      tot_n_points=tot_n_points+1  		! south pole point
      
      ALLOCATE( point_list(1:tot_n_points) )

! loop over theta and phi creating the points required for the triangulation

! north pole point      
      point_count=1
      theta=0d0
      phi=0d0
      
      radius=radius0
      if (problem_surfaces(surface_number)%roughness_flag) then
        CALL random_number(r_random)
	r_random=(r_random-0.5d0)*2d0
        radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
      end if
      CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      point_list(point_count)=point1
     
      do theta_loop=2,n_theta
      
        do phi_loop=1,n_phi
	
	  point_count=point_count+1
	
	  theta=(theta_loop-1)*dtheta	  
	  phi=(phi_loop-1)*dphi
      
          radius=radius0
          if (problem_surfaces(surface_number)%roughness_flag) then
            CALL random_number(r_random)
	    r_random=(r_random-0.5d0)*2d0
            radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
          end if
          CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
          CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
          point_list(point_count)=point1
	
	end do ! next phi
	
      end do ! next theta
	
! south pole point      
      point_count=point_count+1

      theta=pi
      phi=0d0
     
      radius=radius0
      if (problem_surfaces(surface_number)%roughness_flag) then
     	CALL random_number(r_random)
	r_random=(r_random-0.5d0)*2d0
     	radius=radius+r_random*problem_surfaces(surface_number)%roughness_p1
      end if
      CALL rthetaphi_to_xyz_point(radius,theta,phi,point1)
      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      point_list(point_count)=point1

      if (point_count.NE.tot_n_points) then
        write(*,*)'Error in build_surface_sphere'
	write(*,*)point_count,tot_n_points
	STOP
      end if
      
! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=n_theta*n_phi*2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! loop over theta and phi creating surface triangles
      
      triangle_count=0
      
      do theta_loop=1,n_theta
      
        do phi_loop=1,n_phi
	
	  CALL get_point_from_list(theta_loop  ,phi_loop  ,n_theta,n_phi,point1,point_list,tot_n_points)
	  CALL get_point_from_list(theta_loop  ,phi_loop+1,n_theta,n_phi,point2,point_list,tot_n_points)
	  CALL get_point_from_list(theta_loop+1,phi_loop+1,n_theta,n_phi,point3,point_list,tot_n_points)
	  CALL get_point_from_list(theta_loop+1,phi_loop  ,n_theta,n_phi,point4,point_list,tot_n_points)

! set triangles to make the normal point outwards	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point3
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point2
	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point4
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3
	
	end do ! next phi
	
      end do ! next theta
      
      DEALLOCATE( point_list )

  RETURN

END SUBROUTINE build_surface_sphere_OLD

!
!
! NAME
!     SUBROUTINE build_surface_sphere_OLD_OLD
!
! DESCRIPTION
!     build_surface_sphere:
!
!     create a triangulated sphere surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_surface_sphere_OLD_OLD(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: theta1,phi1
real*8	:: theta2,phi2
real*8	:: dtheta,dphi
real*8	:: radius
real*8	:: circumference
integer	:: n_circumference
integer :: n_theta,n_phi
integer :: theta_loop,phi_loop

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3,point4
type(xyz)	:: point5,point6,point7,point8

! START

    
! CREATE A TRIANGULATED SPHERE MESH ON THE SCALE OF DL.    
      
      radius=problem_surfaces(surface_number)%surface_parameters(1)

! set the length of triangle edges to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 12 edges around the circumference

      circumference=2d0*pi*radius
      n_circumference=2*NINT(circumference/dl)  ! ensure that n_circumference is an even number
      
      if (n_circumference.lt.12) n_circumference=12
      
      n_theta=n_circumference/2
      n_phi=n_circumference
      
      dtheta=pi/n_theta
      dphi=2d0*pi/n_phi

! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=n_theta*n_phi*2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! loop over theta and phi creating surface triangles
      
      triangle_count=0
      
      do theta_loop=1,n_theta
      
        do phi_loop=1,n_phi
	
	  theta1=(theta_loop-1)*dtheta
	  theta2= theta_loop   *dtheta
	  
	  phi1=(phi_loop-1)*dphi
	  phi2= phi_loop   *dphi
	  if (phi_loop.eq.1) phi1=0d0
	  if (phi_loop.eq.n_phi) phi2=0d0

! get four points at the corners of a quad patch on the sphere surface
	  
	  CALL rthetaphi_to_xyz_point(radius,theta1,phi1,point1)
	  CALL rthetaphi_to_xyz_point(radius,theta1,phi2,point2)
	  CALL rthetaphi_to_xyz_point(radius,theta2,phi2,point3)
	  CALL rthetaphi_to_xyz_point(radius,theta2,phi1,point4)

! apply the transformation to each of the surface points
          CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
          CALL apply_transformation(point2,problem_surfaces(surface_number)%trans)
          CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
          CALL apply_transformation(point4,problem_surfaces(surface_number)%trans)

! set triangles to make the normal point outwards	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point3
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point2
	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point4
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3
	
	end do ! next phi
	
      end do ! next theta

  RETURN

END SUBROUTINE build_surface_sphere_OLD_OLD
