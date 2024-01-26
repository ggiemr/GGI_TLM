!SUBROUTINE build_surface_stl_triangulated_surface
!SUBROUTINE read_vertex_coordinates!
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
!SUBROUTINE build_surface_stl_triangulated_surface
!
! NAME
!     SUBROUTINE build_surface_stl_triangulated_surface
!
! DESCRIPTION
!     build_surface_stl_triangulated_surface:
!
!     read triangulated_surface from stl format file
!
! Format is:
!  solid name
!  facet normal 5.972882e-015 -3.877537e-002 9.992480e-001
!    outer loop
!      vertex   -2.692400e+000 4.231861e-001 -5.444578e+000
!      vertex   -2.840990e+000 0.000000e+000 -5.461000e+000
!      vertex   -2.692400e+000 0.000000e+000 -5.461000e+000
!    endloop
!  endfacet
!  facet normal 5.987367e-015 -3.877537e-002 9.992480e-001
!    outer loop
!      vertex   -2.692400e+000 4.231861e-001 -5.444578e+000
!      vertex   -2.840990e+000 4.231861e-001 -5.444578e+000
!      vertex   -2.840990e+000 0.000000e+000 -5.461000e+000
!    endloop
!  endfacet
!  .
!  .
!  .
!endsolid
!
! COMMENTS
!  Do we have to read and check the triangle normal direction or can we rely on the ordering of the 
!  vertives to specify this? At the moment we assume that ther vertex order specifies the normal     
!
! HISTORY
!
!     started 2/12/2014 CJS based on build_surface_vtk_triangulated_surface
!     7/1/2019 CJS correct reversal of surface normal by changing the triangle vertex order
!     4/9/2019 CJS correct input format error to generalise things
!
SUBROUTINE build_surface_stl_triangulated_surface(surface_number)

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

integer	:: line_number
integer :: found

real	:: x,y,z

integer	:: n_triangles

integer	:: point,triangle

type(xyz)	:: triangle_point1,triangle_point2,triangle_point3

character*256 	:: ipline

! START

  scale_factor=problem_surfaces(surface_number)%surface_parameters(1)
  reverse_normal=problem_surfaces(surface_number)%surface_parameters(2)

  OPEN(UNIT=local_file_unit,file= problem_surfaces(surface_number)%filename)

! Read the file to check the format is as expected and count the number of triangles 
  n_triangles=0

  line_number=0
  
  read(local_file_unit,'(A)',err=9000)ipline
  line_number=line_number+1
  call convert_to_lower_case(ipline,256)
  found=INDEX(ipline,'solid')
  if (found.eq.0) then
! Format error...
    GOTO 9010
  end if 
  
! Triangle reading loop
10  CONTINUE
    read(local_file_unit,'(A)',err=9000)ipline
    line_number=line_number+1
    call convert_to_lower_case(ipline,256)

! check for the end of the triangle list    
    found=INDEX(ipline,'endsolid')
    if (found.NE.0) then
! end of triangle list
      GOTO 100
    end if

! check for a valid triangle        
    found=INDEX(ipline,'facet')
    if (found.NE.0) then
! read a triangle

      read(local_file_unit,'(A)',err=9000)ipline
      line_number=line_number+1
      call convert_to_lower_case(ipline,256)
      
      found=INDEX(ipline,'outer loop')
      if (found.NE.0) then
! proceed to read three vertices
 
        do point=1,3

          read(local_file_unit,'(A)',err=9000)ipline
          line_number=line_number+1
          call convert_to_lower_case(ipline,256)
	  CALL read_vertex_coordinates(ipline,'vertex',x,y,z,found)

        end do
      
      else
! something unexpected in the file...

        GOTO 9010 
	
      end if ! reading triangle
      
! finish reading triangle
      read(local_file_unit,'(A)',err=9000)ipline
      line_number=line_number+1
      call convert_to_lower_case(ipline,256)
      found=INDEX(ipline,'endloop')
      if (found.EQ.0) then
! Format error...
        GOTO 9010
      end if
      
! finish reading triangle
      read(local_file_unit,'(A)',err=9000)ipline
      line_number=line_number+1
      call convert_to_lower_case(ipline,256)
      found=INDEX(ipline,'endfacet')
      if (found.EQ.0) then
! Format error...
        GOTO 9010
      end if

! Finished reading triangle so can count it            
      n_triangles=n_triangles+1

    else
! something unexpected in the file...

      GOTO 9010 
    
    end if 

  GOTO 10 ! Read next triangle

100 CONTINUE
    
! We have read the file successfully so go bak and read the triangle vertex data
    
  rewind(unit=local_file_unit)
    
  problem_surfaces(surface_number)%number_of_triangles=n_triangles   
  allocate( problem_surfaces(surface_number)%triangle_list(1:n_triangles ) )

! we have checked the file format so we can read it without further checks here

  read(local_file_unit,'(A)',err=9000)ipline  
  
  do triangle=1,n_triangles

! read two lines
    read(local_file_unit,'(A)',err=9000)ipline  
    read(local_file_unit,'(A)',err=9000)ipline  

! read three vertices
    
    read(local_file_unit,'(A)',err=9000)ipline
    line_number=line_number+1
    call convert_to_lower_case(ipline,256)
    CALL read_vertex_coordinates(ipline,'vertex',x,y,z,found)
    triangle_point1%x=x*scale_factor
    triangle_point1%y=y*scale_factor
    triangle_point1%z=z*scale_factor
    
    read(local_file_unit,'(A)',err=9000)ipline
    line_number=line_number+1
    call convert_to_lower_case(ipline,256)
    CALL read_vertex_coordinates(ipline,'vertex',x,y,z,found)
    triangle_point2%x=x*scale_factor
    triangle_point2%y=y*scale_factor
    triangle_point2%z=z*scale_factor
    
    read(local_file_unit,'(A)',err=9000)ipline
    line_number=line_number+1
    call convert_to_lower_case(ipline,256)
    CALL read_vertex_coordinates(ipline,'vertex',x,y,z,found)
    triangle_point3%x=x*scale_factor
    triangle_point3%y=y*scale_factor
    triangle_point3%z=z*scale_factor
    
! apply the transformation to each of the points
    CALL apply_transformation(triangle_point1,problem_surfaces(surface_number)%trans)
    CALL apply_transformation(triangle_point2,problem_surfaces(surface_number)%trans)
    CALL apply_transformation(triangle_point3,problem_surfaces(surface_number)%trans)
        
    if (reverse_normal.ge.0d0) then    ! maintain orientation
      problem_surfaces(surface_number)%triangle_list(triangle)%vertex(1)=triangle_point1
      problem_surfaces(surface_number)%triangle_list(triangle)%vertex(2)=triangle_point2
      problem_surfaces(surface_number)%triangle_list(triangle)%vertex(3)=triangle_point3
    else                               ! reverse orientation
      problem_surfaces(surface_number)%triangle_list(triangle)%vertex(1)=triangle_point1
      problem_surfaces(surface_number)%triangle_list(triangle)%vertex(2)=triangle_point3
      problem_surfaces(surface_number)%triangle_list(triangle)%vertex(3)=triangle_point2 
    end if

! read two lines
    read(local_file_unit,'(A)',err=9000)ipline  
    read(local_file_unit,'(A)',err=9000)ipline  

  end do ! next triangle
  
  CLOSE(UNIT=local_file_unit)
  
  RETURN

9000 CALL write_line('Error in build_surface_stl_triangulated_surface:',0,.TRUE.)
     CALL write_line_integer('Error reading line number:',line_number,0,.TRUE.)
     STOP

9010 CALL write_line('Error in build_surface_stl_triangulated_surface:',0,.TRUE.)
     CALL write_line_integer('Format error line number:',line_number,0,.TRUE.)
     STOP


END SUBROUTINE build_surface_stl_triangulated_surface
!
! ______________________________________________________________________________
!
!
SUBROUTINE read_vertex_coordinates(ipline,test_string,x,y,z,found)

IMPLICIT NONE

character*256	:: ipline
character(*) 	:: test_string

real	:: x,y,z
integer :: found

! local variables

character*256 	:: stripped_ipline
integer first_char,last_char

! START

found=INDEX(ipline,trim(test_string))

if (found.EQ.0) then
  x=0d0
  y=0d0
  z=0d0
  RETURN
end if

first_char=found+len(trim(test_string))
last_char=len(trim(ipline))

stripped_ipline=ipline(first_char:last_char)

read(stripped_ipline,*)x,y,z

RETURN

END SUBROUTINE read_vertex_coordinates
