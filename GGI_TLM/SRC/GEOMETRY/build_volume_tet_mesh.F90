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
!SUBROUTINE build_volume_tet_mesh
!
! NAME
!     SUBROUTINE build_volume_tet_mesh
!
! DESCRIPTION
!     build_volume_tet_mesh:
!
!     read tet mesh from file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 29/04/2013 CJS
!
!
SUBROUTINE build_volume_tet_mesh(volume_number)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: volume_number

! local variables

real*8	:: scale_factor

integer	:: n_tet_points
integer	:: max_point_number

integer	:: n_tets
integer	:: max_tet_number

integer	:: point,tet    
integer	:: point1,point2,point3,point4
integer	:: ip_point
integer	:: ip_tet

real*8	:: x,y,z

integer	:: n_OK_tets
integer	:: tet_count

type(xyz)	:: tet_point1,tet_point2,tet_point3,tet_point4

real*8,allocatable	:: points(:,:)
integer,allocatable	:: tets(:,:)
logical,allocatable	:: tet_OK(:)

character*256 	:: ipline

! START

  scale_factor=problem_volumes(volume_number)%volume_parameters(1)

  OPEN(UNIT=local_file_unit,file= problem_volumes(volume_number)%filename)

! work out the number of points  
  n_tet_points=0
  max_point_number=0
  n_tets=0
  max_tet_number=0

! read up to the first coordinate  
10  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:11).NE.'coordinates') then
      goto 10
    end if 

! read coordinates    

    CALL write_line('starting to read coordinates',0,.TRUE.)
20  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:15).NE.'end coordinates') then
      n_tet_points=n_tet_points+1
      goto 20
    end if
    
  CALL write_line_integer('Number of points=',n_tet_points,0,.TRUE.)
    
! read up to the first tet  
30  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:8).NE.'elements') then
      goto 30
    end if 

! read tets 
40  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:12).NE.'end elements') then
      n_tets=n_tets+1
      goto 40
    end if
    
  CALL write_line_integer('Number of tets=',n_tets,0,.TRUE.)
    
   rewind (local_file_unit)
  
! read the maximum value in the point numbering
! read up to the first coordinate  
50  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:11).NE.'coordinates') then
      goto 50
    end if 

  do point=1,n_tet_points
    read(local_file_unit,*,err=9000)ip_point
    if (ip_point.gt.max_point_number) max_point_number=ip_point
  end do
    
  CALL write_line_integer('Maximum point number=',max_point_number,0,.TRUE.)
  
! read up to the first tet  
55  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:8).NE.'elements') then
      goto 55
    end if 
    
  do tet=1,n_tets
    read(local_file_unit,*,err=9000)ip_tet
    if (ip_tet.gt.max_tet_number) max_tet_number=ip_tet
  end do
  
  CALL write_line_integer('Maximum tet number=',max_tet_number,0,.TRUE.)
  
  rewind (local_file_unit)

! allocate memory for geometry data 
    
  CALL write_line('Allocating memory',0,.TRUE.)
  
  ALLOCATE ( points(1:max_point_number,1:3) ) 
  ALLOCATE ( tets(1:max_tet_number,1:4) ) 
  ALLOCATE ( tet_OK(1:max_tet_number) ) 
  
  points(1:max_point_number,1:3)=0d0
  tets(1:max_tet_number,1:4)=0
  tet_OK(1:max_tet_number)=.FALSE.
  
! read up to the first coordinate  
60  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:11).NE.'coordinates') then
      goto 60
    end if 

! read coordinates
  CALL write_line('Reading points',0,.TRUE.)

  do point=1,n_tet_points
  
    read(local_file_unit,*,err=9000)ip_point,x,y,z
    
    points(ip_point,1)=x*scale_factor
    points(ip_point,2)=y*scale_factor
    points(ip_point,3)=z*scale_factor
  end do
    
! read up to the first tet  
70  CONTINUE

    read(local_file_unit,'(A)',err=9000)ipline
    call convert_to_lower_case(ipline,256)
    if (ipline(1:8).NE.'elements') then
      goto 70
    end if 
    
! read tets    
  CALL write_line('Reading elements',0,.TRUE.)

  do tet=1,n_tets
  
    read(local_file_unit,*,err=9000)ip_tet,point1,point2,point3,point4
    
    tets(ip_tet,1)=point1
    tets(ip_tet,2)=point2
    tets(ip_tet,3)=point3
    tets(ip_tet,4)=point4
    
  end do
  
  CLOSE(UNIT=local_file_unit)

! Check the number of OK tets

  n_OK_tets=0
  
  do tet=1,n_tets

    point1=tets(tet,1)
    point2=tets(tet,2)
    point3=tets(tet,3)
    point4=tets(tet,4)
    
    if ( (point1.ne.0).AND.(point2.ne.0).AND.(point3.ne.0).AND.(point4.ne.0) ) then
! we have a valid tet so set the cell data for this tet

      n_OK_tets=n_OK_tets+1
      
    end if

  end do
  
  problem_volumes(volume_number)%number_of_tets=n_OK_tets   
  allocate( problem_volumes(volume_number)%tet_list(1:n_OK_tets ) )
  
  tet_count=0
  
  do tet=1,n_tets

    point1=tets(tet,1)
    point2=tets(tet,2)
    point3=tets(tet,3)
    point4=tets(tet,4)
    
    if ( (point1.ne.0).AND.(point2.ne.0).AND.(point3.ne.0).AND.(point4.ne.0) ) then
! we have a valid tet so set the cell data for this tet
      
      tet_point1%x=points(point1,1)
      tet_point1%y=points(point1,2)
      tet_point1%z=points(point1,3)
      tet_point2%x=points(point2,1)
      tet_point2%y=points(point2,2)
      tet_point2%z=points(point2,3)
      tet_point3%x=points(point3,1)
      tet_point3%y=points(point3,2)
      tet_point3%z=points(point3,3)
      tet_point4%x=points(point4,1)
      tet_point4%y=points(point4,2)
      tet_point4%z=points(point4,3)
    
! apply the transformation to each of the points
      CALL apply_transformation(tet_point1,problem_volumes(volume_number)%trans)
      CALL apply_transformation(tet_point2,problem_volumes(volume_number)%trans)
      CALL apply_transformation(tet_point3,problem_volumes(volume_number)%trans)
      CALL apply_transformation(tet_point4,problem_volumes(volume_number)%trans)
	  
      tet_count=tet_count+1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(1)=tet_point1
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(2)=tet_point2
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(3)=tet_point3
      problem_volumes(volume_number)%tet_list(tet_count)%vertex(4)=tet_point4
           
    end if
      

  end do
   
  if (allocated(points)) DEALLOCATE ( points ) 
  if (allocated(tets))   DEALLOCATE ( tets ) 
  if (allocated(tet_OK)) DEALLOCATE ( tet_OK ) 
  

  RETURN

9000 CALL write_line('Error in build_volume_tet_mesh:',0,.TRUE.)
     CALL write_line_integer('volume_number=',volume_number,0,.TRUE.)
     STOP


END SUBROUTINE build_volume_tet_mesh
