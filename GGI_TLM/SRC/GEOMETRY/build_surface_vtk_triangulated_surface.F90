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
!SUBROUTINE build_surface_vtk_triangulated_surface
!
! NAME
!     SUBROUTINE build_surface_vtk_triangulated_surface
!
! DESCRIPTION
!     build_surface_vtk_triangulated_surface:
!
!     read triangulated_surface from vtk format file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 22/04/2013 CJS based on build_surface_triangulated_surface
!
!
SUBROUTINE build_surface_vtk_triangulated_surface(surface_number)

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

integer	:: n_triangles

integer :: count,triangle_count
integer :: point1,point2,point3

integer	:: point,triangle

type(xyz)	:: triangle_point1,triangle_point2,triangle_point3

real*8,allocatable	:: input_reals(:)
real*8,allocatable	:: points(:,:)
integer,allocatable	:: triangles(:,:)

character*256 	:: ipline

! START

  scale_factor=problem_surfaces(surface_number)%surface_parameters(1)
  reverse_normal=problem_surfaces(surface_number)%surface_parameters(2)

  OPEN(UNIT=local_file_unit,file= problem_surfaces(surface_number)%filename)

! work out the number of points  
  n_triangle_points=0
  n_triangles=0

! read up to the first coordinate  
10  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:6).NE.'points') then
      GOTO 10
    end if 
    
  backspace(unit=local_file_unit)
    
  read(local_file_unit,*)ipline(1:6),n_triangle_points
    
  CALL write_line_integer('Number of points=',n_triangle_points,0,.TRUE.)

! read coordinates    
  
  ALLOCATE ( input_reals(1:n_triangle_points*3) ) 
  
  read(local_file_unit,*)input_reals(1:n_triangle_points*3)
  
  ALLOCATE ( points(1:n_triangle_points,1:3) ) 
  
  count=0
  
  do point=1,n_triangle_points
  
    count=count+1
    points(point,1)=input_reals(count)*scale_factor
    count=count+1
    points(point,2)=input_reals(count)*scale_factor
    count=count+1
    points(point,3)=input_reals(count)*scale_factor
    
  end do

  CALL write_line('Finished reading points=',0,.TRUE.)

    CALL write_line('Starting to read polygons',0,.TRUE.)
20  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:8).NE.'polygons') then
      GOTO 20
    end if
       
  backspace(unit=local_file_unit)
    
  read(local_file_unit,*)ipline(1:8),n_triangles
    
  CALL write_line_integer('Number of polygons=',n_triangles,0,.TRUE.)
    
! read triangles    
  CALL write_line('Reading triangles',0,.TRUE.)
  
  ALLOCATE ( triangles(1:n_triangles,1:3) ) 

  do triangle=1,n_triangles
  
    read(local_file_unit,*,err=9000)count,point1,point2,point3
    
    if (count.ne.3) GOTO 9010
    
    if (reverse_normal.ge.0d0) then ! maintain orientation
      triangles(triangle,1)=point1+1
      triangles(triangle,2)=point2+1
      triangles(triangle,3)=point3+1
    else        				! reverse orientation
      triangles(triangle,1)=point3+1
      triangles(triangle,2)=point2+1
      triangles(triangle,3)=point1+1  
    end if
  end do
  
  CLOSE(UNIT=local_file_unit)

! Transfer data onto the correct surface structure
  
  problem_surfaces(surface_number)%number_of_triangles=n_triangles   
  allocate( problem_surfaces(surface_number)%triangle_list(1:n_triangles ) )
  
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
   
  if (allocated(input_reals)) DEALLOCATE ( input_reals ) 
  if (allocated(points)) DEALLOCATE ( points ) 
  if (allocated(triangles)) DEALLOCATE ( triangles ) 
 

  RETURN

9000 CALL write_line('Error in build_surface_vtk_triangulated_surface:',0,.TRUE.)
     CALL write_line_integer('surface_number=',surface_number,0,.TRUE.)
     STOP

9010 CALL write_line('Error in build_surface_vtk_triangulated_surface:',0,.TRUE.)
     CALL write_line_integer('Polygon is not a triangle: ',triangle,0,.TRUE.)
     STOP


END SUBROUTINE build_surface_vtk_triangulated_surface
