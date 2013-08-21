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
!SUBROUTINE read_mesh
!
! NAME
!     SUBROUTINE read_mesh
!
! DESCRIPTION
!     read_mesh:
!
!     read volume, surface, line and point meshes from file
!
!     The parallel implementation only reads the volume, line and surface information
!     within the current process. All points are read into all processes as these are
!     defined on a slightly different basis - this maybe something to think about changing...
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 12/08/2012 CJS
!     parallel 22/11/2012
!
SUBROUTINE read_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants
USE cell_parameters

IMPLICIT NONE

! local variables

integer	:: volume_number
integer	:: cell
integer :: total_number_of_cells
integer :: cell_count

integer	:: surface_number
integer	:: face
integer :: total_number_of_faces
integer :: face_count

integer	:: line_number
integer	:: segment
integer :: total_number_of_cell_segments
integer :: cell_segment_count

integer	:: point_number
integer	:: point
integer	:: total_number_of_points
integer	:: point_count

integer :: read_loop

integer :: i1,j1,k1,point1
integer :: i2,j2,k2,point2

real*8	:: x,y,z

! START

  CALL write_line('CALLED: read_mesh',0,output_to_screen_flag)

! Open mesh file

  CALL open_file(mesh_file_unit,mesh_file_extension)

  do read_loop=1,2
  
    rewind(unit=mesh_file_unit)

! read general mesh parameters

    read(mesh_file_unit,*)dl
    read(mesh_file_unit,*)nx,ny,nz
    read(mesh_file_unit,*)mesh_xmin,mesh_xmax
    read(mesh_file_unit,*)mesh_ymin,mesh_ymax
    read(mesh_file_unit,*)mesh_zmin,mesh_zmax
  
    if (read_loop.eq.1) then
! Work out the size of the mesh required for this processor
      CALL mesh_partition()
    end if

! STAGE 1: READ VOLUMES
    
    read(mesh_file_unit,*)n_volumes
  
    if (read_loop.eq.1) then 
! Allocate the structure required for each volume
      ALLOCATE( problem_volumes(1:n_volumes) )
    end if
    
! read the cell list for each of the volumes
    
    do volume_number=1,n_volumes
    
      read(mesh_file_unit,*)total_number_of_cells
      
      cell_count=0
          
      do cell=1,total_number_of_cells
! read cell data      
        read(mesh_file_unit,*)	i1,j1,k1,point1

! check whether the cell belongs to this process	
	if ( rank.EQ.cell_rank(k1) ) then
	
	  cell_count=cell_count+1
	  
	  if (read_loop.eq.2) then
! put the data into the already allocated cell_list	  
	    problem_volumes(volume_number)%cell_list(cell_count)%cell%i=i1
	    problem_volumes(volume_number)%cell_list(cell_count)%cell%j=j1
	    problem_volumes(volume_number)%cell_list(cell_count)%cell%k=k1
	    problem_volumes(volume_number)%cell_list(cell_count)%point =point1
	  end if ! read_loop.eq.2
	  
	end if ! cell belongs to this process
	
      end do ! next cell
      
      if (read_loop.eq.1) then
! allocate the cell_list according to the number of cells in this processors mesh

        problem_volumes(volume_number)%number_of_cells=cell_count   
        if (cell_count.gt.0) then
          ALLOCATE( problem_volumes(volume_number)%cell_list(1:cell_count) )
        end if
      
      end if ! read_loop.eq.1
    
    end do ! next volume number

! STAGE 2: READ SURFACES
    
    read(mesh_file_unit,*)n_surfaces
  
    if (read_loop.eq.1) then 
! Allocate the structure required for each surface
      ALLOCATE( problem_surfaces(1:n_surfaces) )
    end if
    
! read the face list for each of the surfaces
    
    do surface_number=1,n_surfaces
    
      read(mesh_file_unit,*)total_number_of_faces

      read(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_xmin,problem_surfaces(surface_number)%mesh_xmax
      read(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_ymin,problem_surfaces(surface_number)%mesh_ymax
      read(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_zmin,problem_surfaces(surface_number)%mesh_zmax

      read(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_xmin,problem_surfaces(surface_number)%mesh_cell_xmax
      read(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_ymin,problem_surfaces(surface_number)%mesh_cell_ymax
      read(mesh_file_unit,*)problem_surfaces(surface_number)%mesh_cell_zmin,problem_surfaces(surface_number)%mesh_cell_zmax
      
      face_count=0
          
      do face=1,total_number_of_faces
! read cell data      
        read(mesh_file_unit,*)	i1,j1,k1,point1

! check whether the cell belongs to this process
	
!	if (    ( rank.EQ.cell_rank(k1) )			&  ! OLD TEST
!	    .OR.( (k1.eq.nzmax).AND.(point1.eq.face_zmin) )  	&
!	    .OR.( (k1.eq.nzmin).AND.(point1.eq.face_zmax) )  	&
!	    ) then
	    
	if ( rank.EQ.cell_face_rank(k1,point1) ) then
	
	  face_count=face_count+1
	  
	  if (read_loop.eq.2) then
! put the data into the already allocated cell_list	  
	    problem_surfaces(surface_number)%face_list(face_count)%cell%i=i1
	    problem_surfaces(surface_number)%face_list(face_count)%cell%j=j1
	    problem_surfaces(surface_number)%face_list(face_count)%cell%k=k1
	    problem_surfaces(surface_number)%face_list(face_count)%point =point1
	  end if ! read_loop.eq.2
	  
	end if ! face belongs to this process
	
      end do ! next face
      
      if (read_loop.eq.1) then
! allocate the cell_list according to the number of cells in this processors mesh

        problem_surfaces(surface_number)%number_of_faces=face_count   
        if (face_count.gt.0) then
	
          ALLOCATE( problem_surfaces(surface_number)%face_list(1:face_count) )
        end if
      
      end if ! read_loop.eq.1
    
    end do ! next surface number

! STAGE 3: READ LINES
 
    read(mesh_file_unit,*)n_lines
  
    if (read_loop.eq.1) then 
! Allocate the structure required for each line
      ALLOCATE( problem_lines(1:n_lines) )
    end if
    
! read the line list for each of the lines
    
    do line_number=1,n_lines
    
      read(mesh_file_unit,*)total_number_of_cell_segments
      
      cell_segment_count=0
          
      do segment=1,total_number_of_cell_segments
      
! read cell_segment data      
        read(mesh_file_unit,*)	i1,j1,k1,point1,i2,j2,k2,point2

! check whether the cell_segment belongs to this process	
	if ( rank.EQ.cell_rank(k1)  ) then
	
	  cell_segment_count=cell_segment_count+1
	  
	  if (read_loop.eq.2) then
! put the data into the already allocated cell_segment_list	  
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(1)%cell%i=i1
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(1)%cell%j=j1
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(1)%cell%k=k1
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(1)%point= point1
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(2)%cell%i=i2
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(2)%cell%j=j2
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(2)%cell%k=k2
	    problem_lines(line_number)%cell_segment_list(cell_segment_count)%segment_point(2)%point= point2
	  end if ! read_loop.eq.2
	  
	end if ! cell_segment belongs to this process
	
      end do ! next cell_segment
      
      if (read_loop.eq.1) then
! allocate the cell_segment_list according to the number of cell_segments in this processors mesh

        problem_lines(line_number)%number_of_cell_segments=cell_segment_count   
        if (cell_segment_count.gt.0) then
          ALLOCATE( problem_lines(line_number)%cell_segment_list(1:cell_segment_count) )
        end if
      
      end if ! read_loop.eq.1
    
    end do ! next line number

! STAGE 4: READ POINTS

    read(mesh_file_unit,*)total_number_of_points
    
 ! read the point list 

    point_count=0          
    do point=1,total_number_of_points
! read point data      

      read(mesh_file_unit,*)	x,y,z
      read(mesh_file_unit,*)	i1,j1,k1
      read(mesh_file_unit,*)	i2,j2,k2,point2

!Note: All points are read into all processes ...	
!**** if ( rank.EQ.cell_rank(k1) ) then
	
	point_count=point_count+1
	  
	if (read_loop.eq.2) then
! put the data into the already allocated point_list	  
	  problem_points(point_count)%point%x=x
	  problem_points(point_count)%point%y=y
	  problem_points(point_count)%point%z=z
	  problem_points(point_count)%cell%i=i1
	  problem_points(point_count)%cell%j=j1
	  problem_points(point_count)%cell%k=k1
	  problem_points(point_count)%face%cell%i=i2
	  problem_points(point_count)%face%cell%j=j2
	  problem_points(point_count)%face%cell%k=k2
	  problem_points(point_count)%face%point =point2
	end if ! read_loop.eq.2
	  
!**** end if ! point belongs to this process
	
    end do ! next point
      
    if (read_loop.eq.1) then
! allocate the point_list according to the number of points in this processors mesh

      n_points=point_count   
      if (n_points.gt.0) then
        ALLOCATE( problem_points(n_points) )
      end if
      
    end if ! read_loop.eq.1
    
  end do ! next read_loop
  
  CALL close_file(mesh_file_unit)

  CALL write_line('FINISHED: read_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE read_mesh
