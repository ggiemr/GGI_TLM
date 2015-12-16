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
!SUBROUTINE build_surface_rectangle
!SUBROUTINE build_surface_rectangle_OLD
!
! NAME
!     SUBROUTINE build_surface_rectangle
!
! DESCRIPTION
!     build_surface_rectangle:
!
!     create a triangulated rectangle surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!     surface roughness 18/8/2015 CJS
!
!
SUBROUTINE build_surface_rectangle(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: lx,ly
      
real*8	:: x,y,z
real*8	:: dx,dy
real*8	:: nlx,nly
integer :: n_x,n_y
integer :: x_loop,y_loop

integer	:: number_of_triangles
integer	:: triangle_count

type(xyz)	:: point1,point2,point3,point4

integer		:: tot_n_points
integer		:: point_count
type(xyz),allocatable	:: point_list(:)
real*8		:: r_random
integer		:: tvalue
integer		:: pvalue

real*8		:: edge_length

! START

    
! CREATE A TRIANGULATED SPHERE MESH ON THE SCALE OF DL.    
          
! create rectangle vertices
      lx=problem_surfaces(surface_number)%surface_parameters(1)
      ly=problem_surfaces(surface_number)%surface_parameters(2)

! set the length of triangle edges based on the correlation length of the surface roughness parameter 
! if this is unset or too small then set it to be of the order of dl/2, half the mesh edge length and ensure
! we have at least 2 edges in each dimension

      edge_length=problem_surfaces(surface_number)%roughness_p2*2d0
      if (edge_length.LT.dl/2d0) edge_length=dl
      
      nlx=2*NINT(lx/edge_length)  ! ensure that nlx is an even number
      if (nlx.lt.2) nlx=2
      nly=2*NINT(ly/edge_length)  ! ensure that nlx is an even number
      if (nly.lt.2) nly=2
      
      n_x=nlx
      n_y=nly
      
      dx=lx/(n_x-1)
      dy=ly/(n_y-1)
      
! calculate the number of points and allocate memory for the point data
      tot_n_points=n_x*n_y  	! surface points
      
      ALLOCATE( point_list(1:tot_n_points) )

! loop over x and y creating the points required for the triangulation

      point_count=0
           
      do x_loop=1,n_x
      
        do y_loop=1,n_y
	
	  point_count=point_count+1
	
	  x=(x_loop-1)*dx-lx/2d0
	  y=(y_loop-1)*dy-ly/2d0
          
	  if ( (x_loop.NE.1).AND.(x_loop.NE.n_x).AND.(y_loop.NE.1).AND.(y_loop.NE.n_y) ) then
! add surface roughness if required
            if (problem_surfaces(surface_number)%roughness_flag) then
              CALL random_number(r_random)
	      r_random=(r_random-0.5d0)*2d0
              z=r_random*problem_surfaces(surface_number)%roughness_p1
            end if
	  else
	    z=0d0
	  end if
	  
          point1%x=x
          point1%y=y
          point1%z=z
          CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
          point_list(point_count)=point1
	
	end do ! next y
	
      end do ! next x
	
      if (point_count.NE.tot_n_points) then
        write(*,*)'Error in build_surface_rectangle'
	write(*,*)point_count,tot_n_points
	STOP
      end if
      
! calculate the number of triangles and allocate memory for the triangulated surface data      
      number_of_triangles=(n_x-1)*(n_y-1)*2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! loop over x and y creating surface triangles
      
      triangle_count=0
      
      do x_loop=1,n_x-1
      
        do y_loop=1,n_y-1
	
          point1=point_list(1+(x_loop-1)*n_y+(y_loop-1))
          point2=point_list(1+(x_loop-1)*n_y+(y_loop  ))
          point3=point_list(1+(x_loop  )*n_y+(y_loop  ))
          point4=point_list(1+(x_loop  )*n_y+(y_loop-1))

! set triangles to make the normal point outwards	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point3
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point2
	  
	  triangle_count=triangle_count+1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point4
	  problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3
	
	end do ! next y
	
      end do ! next x
      
      if (triangle_count.NE.number_of_triangles) then
        write(*,*)'Error in build_surface_rectangle'
	write(*,*)triangle_count,number_of_triangles
	STOP
      end if
      
      DEALLOCATE( point_list )

  RETURN

END SUBROUTINE build_surface_rectangle
!
! NAME
!     SUBROUTINE build_surface_rectangle_OLD
!
! DESCRIPTION
!     build_surface_rectangle:
!
!     create a triangulated rectangle surface
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_surface_rectangle_OLD(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

integer	:: number_of_triangles
integer	:: triangle_count

real*8	:: lx,ly

type(xyz)	:: point1,point2,point3,point4

! START

! CREATE A TRIANGULATED RECTANGULAR BLOCK  
    
      number_of_triangles=2
      problem_surfaces(surface_number)%number_of_triangles=number_of_triangles    
      allocate( problem_surfaces(surface_number)%triangle_list(1:number_of_triangles) )

! create rectangle vertices
      lx=problem_surfaces(surface_number)%surface_parameters(1)
      ly=problem_surfaces(surface_number)%surface_parameters(2)

      point1%x=-lx/2d0
      point1%y=-ly/2d0
      point1%z=0d0
      point2%x=+lx/2d0
      point2%y=-ly/2d0
      point2%z=0d0
      point3%x=+lx/2d0
      point3%y=+ly/2d0
      point3%z=0d0
      point4%x=-lx/2d0
      point4%y=+ly/2d0
      point4%z=0d0
    
! apply the transformation to each of the points
      CALL apply_transformation(point1,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point2,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point3,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(point4,problem_surfaces(surface_number)%trans)
	  
! create surface triangles    
	  
      triangle_count=0
	  
      triangle_count=triangle_count+1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point2
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point3

      triangle_count=triangle_count+1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=point1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=point3
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=point4	       

  RETURN

END SUBROUTINE build_surface_rectangle_OLD
