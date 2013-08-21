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
!SUBROUTINE build_surface_triangulated_surface
!
! NAME
!     SUBROUTINE build_surface_triangulated_surface
!
! DESCRIPTION
!     build_surface_triangulated_surface:
!
!     read triangulated_surface from file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/08/2012 CJS
!
!
SUBROUTINE build_surface_triangulated_surface(surface_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: surface_number

! local variables

real*8	:: scale_factor
real*8	:: reverse_normal

integer	:: n_triangle_points
integer	:: max_point_number

integer	:: n_triangles
integer	:: max_triangle_number

integer	:: point,triangle    
integer	:: point1,point2,point3
integer	:: ip_point
integer	:: ip_triangle

real*8	:: x,y,z

integer	:: n_OK_triangles
integer	:: triangle_count

type(xyz)	:: triangle_point1,triangle_point2,triangle_point3

real,allocatable	:: points(:,:)
integer,allocatable	:: triangles(:,:)
logical,allocatable	:: triangle_OK(:)

character*256 	:: ipline

! START

  scale_factor=problem_surfaces(surface_number)%surface_parameters(1)
  reverse_normal=problem_surfaces(surface_number)%surface_parameters(2)

  OPEN(UNIT=local_file_unit,file= problem_surfaces(surface_number)%filename)

! work out the number of points  
  n_triangle_points=0
  max_point_number=0
  n_triangles=0
  max_triangle_number=0

! read up to the first coordinate  
10  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:11).NE.'coordinates') then
      goto 10
    end if 

! read coordinates    

    CALL write_line('starting to read coordinates',0,.TRUE.)
20  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:15).NE.'end coordinates') then
      n_triangle_points=n_triangle_points+1
      goto 20
    end if
    
  CALL write_line_integer('Number of points=',n_triangle_points,0,.TRUE.)
    
! read up to the first triangle  
30  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:8).NE.'elements') then
      goto 30
    end if 

! read triangles 
40  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:12).NE.'end elements') then
      n_triangles=n_triangles+1
      goto 40
    end if
    
  CALL write_line_integer('Number of triangles=',n_triangles,0,.TRUE.)
    
   rewind (local_file_unit)
  
! read the maximum value in the point numbering
! read up to the first coordinate  
50  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:11).NE.'coordinates') then
      goto 50
    end if 

  do point=1,n_triangle_points
    read(local_file_unit,*,err=9000)ip_point
    if (ip_point.gt.max_point_number) max_point_number=ip_point
  end do
    
  CALL write_line_integer('Maximum point number=',max_point_number,0,.TRUE.)
  
! read up to the first triangle  
55  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:8).NE.'elements') then
      goto 55
    end if 
    
  do triangle=1,n_triangles
    read(local_file_unit,*,err=9000)ip_triangle
    if (ip_triangle.gt.max_triangle_number) max_triangle_number=ip_triangle
  end do
  
  CALL write_line_integer('Maximum triangle number=',max_triangle_number,0,.TRUE.)
  
  rewind (local_file_unit)

! allocate memory for geometry data 
    
  CALL write_line('Allocating memory',0,.TRUE.)
  
  ALLOCATE ( points(1:max_point_number,1:3) ) 
  ALLOCATE ( triangles(1:max_triangle_number,1:3) ) 
  ALLOCATE ( triangle_OK(1:max_triangle_number) ) 
  
  points(1:max_point_number,1:3)=0d0
  triangles(1:max_triangle_number,1:3)=0
  triangle_OK(1:max_triangle_number)=.FALSE.
  
! read up to the first coordinate  
60  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:11).NE.'coordinates') then
      goto 60
    end if 

! read coordinates
  CALL write_line('Reading points',0,.TRUE.)

  do point=1,n_triangle_points
  
    read(local_file_unit,*,err=9000)ip_point,x,y,z
    
    points(ip_point,1)=x*scale_factor
    points(ip_point,2)=y*scale_factor
    points(ip_point,3)=z*scale_factor
  end do
    
! read up to the first triangle  
70  CONTINUE

    read(local_file_unit,'(A)',err=9000),ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:8).NE.'elements') then
      goto 70
    end if 
    
! read triangles    
  CALL write_line('Reading elements',0,.TRUE.)

  do triangle=1,n_triangles
  
    read(local_file_unit,*,err=9000)ip_triangle,point1,point2,point3
    
    if (reverse_normal.ge.0d0) then ! maintain orientation
      triangles(ip_triangle,1)=point1
      triangles(ip_triangle,2)=point2
      triangles(ip_triangle,3)=point3
    else                                           ! reverse orientation
      triangles(ip_triangle,1)=point3
      triangles(ip_triangle,2)=point2
      triangles(ip_triangle,3)=point1   
    end if
  end do
  
  CLOSE(UNIT=local_file_unit)

! Check the number of OK triangles

  n_OK_triangles=0
  
  do triangle=1,n_triangles

    point1=triangles(triangle,1)
    point2=triangles(triangle,2)
    point3=triangles(triangle,3)
    
    if ( (point1.ne.0).AND.(point2.ne.0).AND.(point3.ne.0) ) then
! we have a valid triangle so set the cell data for this triangle

      n_OK_triangles=n_OK_triangles+1
      
    end if

  end do
  
  problem_surfaces(surface_number)%number_of_triangles=n_OK_triangles   
  allocate( problem_surfaces(surface_number)%triangle_list(1:n_OK_triangles ) )
  
  triangle_count=0
  
  do triangle=1,n_triangles

    point1=triangles(triangle,1)
    point2=triangles(triangle,2)
    point3=triangles(triangle,3)
    
    if ( (point1.ne.0).AND.(point2.ne.0).AND.(point3.ne.0) ) then
! we have a valid triangle so set the cell data for this triangle
      
      triangle_point1%x=points(point1,1)
      triangle_point1%y=points(point1,2)
      triangle_point1%z=points(point1,3)
      triangle_point2%x=points(point2,1)
      triangle_point2%y=points(point2,2)
      triangle_point2%z=points(point2,3)
      triangle_point3%x=points(point3,1)
      triangle_point3%y=points(point3,2)
      triangle_point3%z=points(point3,3)
    
! apply the transformation to each of the points
      CALL apply_transformation(triangle_point1,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(triangle_point2,problem_surfaces(surface_number)%trans)
      CALL apply_transformation(triangle_point3,problem_surfaces(surface_number)%trans)
	  
      triangle_count=triangle_count+1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(1)=triangle_point1
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(2)=triangle_point3
      problem_surfaces(surface_number)%triangle_list(triangle_count)%vertex(3)=triangle_point2
            
    end if
      

  end do
   
  if (allocated(points)) DEALLOCATE ( points ) 
  if (allocated(triangles)) DEALLOCATE ( triangles ) 
  if (allocated(triangle_OK)) DEALLOCATE ( triangle_OK ) 
 

  RETURN

9000 CALL write_line('Error in build_surface_triangulated_surface:',0,.TRUE.)
     CALL write_line_integer('surface_number=',surface_number,0,.TRUE.)
     STOP


END SUBROUTINE build_surface_triangulated_surface
