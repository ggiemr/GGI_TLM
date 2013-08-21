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
!SUBROUTINE create_cable_meshes
!SUBROUTINE filter_cable_segment_list
!
! NAME
!     SUBROUTINE create_cable_meshes
!
! DESCRIPTION
!      Create a cable based mesh list, including the junction points at the
!      cable end points
!      
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!     include terminations to surfaces 19/11/2012 CJS
!
!
SUBROUTINE create_cable_meshes()

USE TLM_general
USE Cables
USE Geometry_types
USE Geometry
USE cell_parameters
USE constants
USE File_information

IMPLICIT NONE

! local variables

  integer	:: cable
  integer	:: junction1,junction2
  integer	:: cable_end_point1,cable_end_point2
  type(ijk)	:: cable_end_cell1,cable_end_cell2
  type(cell_point)	:: cable_end_cell_point1,cable_end_cell_point2
  
  integer		:: line_end_point1,line_end_point2
  type(cell_point)	:: line_end_cell_point1,line_end_cell_point2
  type(cell_point)	:: last_cell_point
  
  type(cell_point)	:: cell_point1,cell_point2
  
  integer	:: line,line_loop,line_number
  
  integer	:: n_line_cell_segments
  integer	:: segment_count
  
  integer			:: local_number_of_cable_segments
  type(cell_segment),allocatable:: local_cable_segment_list(:)
  integer,allocatable		:: local_direction_sign_list(:)
  
  integer	:: cell_segment_number
  integer	:: first_cell_segment
  integer	:: last_cell_segment
  integer	:: cell_segment_step
  integer	:: first_cell_segment_point
  integer	:: second_cell_segment_point
  
  integer	:: segment_sign
  
  integer 	:: i
  logical	:: face_termination
  integer	:: cable_face

! function variables  
  logical	:: same_cell_point
  
! START

  CALL write_line('CALLED: create_cable_meshes',0,output_to_screen_flag)
  
! loop over cables
  do cable=1,n_cables  
  
    write(cable_info_file_unit,*)'Cable',cable
  
! STAGE 1: get end point 1 junction coordinates and cell_points

    junction1=cable_list(cable)%junction_1
    
    if ( (junction1.LT.1).OR.(junction1.GT.n_cable_junctions) ) GOTO 9020
    cable_end_point1=cable_junction_list(junction1)%point_number
    
    if (cable_junction_list(junction1)%junction_type.EQ.junction_type_cell) then
    
      cable_end_cell_point1%cell =problem_points(cable_end_point1)%cell
      cable_end_cell_point1%point=centre

    else
! face junction

! check whether this junction terminaties to a surface and if so ensure that the 
! junction point is on a meshed face. 
      face_termination=.FALSE.
      do i=1,cable_junction_list(junction1)%n_internal_connection_nodes
        if (cable_junction_list(junction1)%BC(i).ne.0) face_termination=.TRUE.
      end do
    
      if (face_termination) then    
        write(cable_info_file_unit,*)'Looking for termination surface, end 1 of cable ',cable
        CALL get_closest_mesh_face(cable_end_point1,cable_face)      
      else ! use the geometrically closest face of the cell    
        cable_face=problem_points(cable_end_point1)%face%point      
      end if

      cable_end_cell_point1%cell =problem_points(cable_end_point1)%cell
      cable_end_cell_point1%point=cable_face

    end if
    
    cable_junction_list(junction1)%cell_point=cable_end_cell_point1
    
    write(cable_info_file_unit,*)'Junction  1',junction1
    write(cable_info_file_unit,*)'end point 1',cable_end_point1
    write(cable_info_file_unit,*)'end cell  1',cable_end_cell_point1
  
! STAGE 2: get end point 2 junction coordinates and cell_points
    
    junction2=cable_list(cable)%junction_2
    
    if ( (junction2.LT.1).OR.(junction2.GT.n_cable_junctions) ) GOTO 9030

    cable_end_point2=cable_junction_list(junction2)%point_number
   
    if (cable_junction_list(junction2)%junction_type.EQ.junction_type_cell) then
    
      cable_end_cell_point2%cell =problem_points(cable_end_point2)%cell
      cable_end_cell_point2%point=centre
      
    else
! face junction

! check whether this junction terminaties to a surface and if so ensure that the 
! junction point is on a meshed face. 
      face_termination=.FALSE.
      do i=1,cable_junction_list(junction2)%n_internal_connection_nodes
        if (cable_junction_list(junction2)%BC(i).ne.0) face_termination=.TRUE.
      end do
    
      if (face_termination) then    
        write(cable_info_file_unit,*)'Looking for termination surface, end 2 of cable ',cable
        CALL get_closest_mesh_face(cable_end_point2,cable_face)      
      else ! use the geometrically closest face of the cell    
        cable_face=problem_points(cable_end_point2)%face%point      
      end if

      cable_end_cell_point2%cell =problem_points(cable_end_point2)%cell
      cable_end_cell_point2%point=cable_face

    end if
    
    cable_junction_list(junction2)%cell_point=cable_end_cell_point2
           
    write(cable_info_file_unit,*)'Junction  2',junction2
    write(cable_info_file_unit,*)'end point 2',cable_end_point2
    write(cable_info_file_unit,*)'end cell  2',cable_end_cell_point2

! STAGE 3: estimate the number of cable segments on this cable route     

    local_number_of_cable_segments=2 ! 2 end point joining segments
    
    do line_loop=1,cable_list(cable)%n_lines   
      line_number=cable_list(cable)%line_list(line_loop)
      local_number_of_cable_segments=local_number_of_cable_segments+problem_lines(line_number)%number_of_cell_segments+1
    end do ! next line_loop (line on the cable route)
    
    write(cable_info_file_unit,*)'cable line list:',cable_list(cable)%line_list(:)
    write(cable_info_file_unit,*)'Estimated number of cell segments',local_number_of_cable_segments
    
! Allocate a local cell segment list for the cable 
    ALLOCATE( local_cable_segment_list(1:local_number_of_cable_segments) )   
     
! Allocate a local direction sign list for the cable 
    ALLOCATE( local_direction_sign_list(1:local_number_of_cable_segments) )    
    
! STAGE 4: set first point of cable to junction 1 cell  

    segment_count=1
    
    local_cable_segment_list(segment_count)%segment_point(1)=cable_end_cell_point1
    
! STAGE 5: if the end point is on a face then create a segment to bring the current point to the cell centre.  

    if (cable_end_cell_point1%point.eq.centre) then
    
      last_cell_point=local_cable_segment_list(segment_count)%segment_point(1)
      segment_count=0
      
    else
    
      local_cable_segment_list(segment_count)%segment_point(2)%cell=	&
            local_cable_segment_list(segment_count)%segment_point(1)%cell
      local_cable_segment_list(segment_count)%segment_point(2)%point=centre
      
      last_cell_point=local_cable_segment_list(segment_count)%segment_point(2)
      local_direction_sign_list(segment_count)=1 ! cable direction is from face to centre
      segment_count=1
	
      write(cable_info_file_unit,*)'Setting segment',segment_count,' direction=',	&
                                   local_direction_sign_list(segment_count)
      write(cable_info_file_unit,*)'      ',local_cable_segment_list(segment_count)%segment_point(1)%cell,   &
                                face_string(local_cable_segment_list(segment_count)%segment_point(1)%point)
      write(cable_info_file_unit,*)'      ',local_cable_segment_list(segment_count)%segment_point(2)%cell,   &
                                face_string(local_cable_segment_list(segment_count)%segment_point(2)%point)
    
    end if
    
! STAGE 6: loop over cable line list    

    segment_sign=-1  ! initial cable direction is from centre to face
    
    do line_loop=1,cable_list(cable)%n_lines 
    
      write(cable_info_file_unit,*)
    
      line=cable_list(cable)%line_list(line_loop)
      
      write(cable_info_file_unit,*)'Line on cable route',line_loop,' line number',line
      
      n_line_cell_segments=problem_lines(line)%number_of_cell_segments
      
      write(cable_info_file_unit,*)'number of cell segments',n_line_cell_segments
      
      if (n_line_cell_segments.gt.0) then
      
! get the line end points
        line_end_point1=1
        line_end_cell_point1=problem_lines(line)%cell_segment_list(line_end_point1)%segment_point(1)
      
        line_end_point2=n_line_cell_segments
        line_end_cell_point2=problem_lines(line)%cell_segment_list(line_end_point2)%segment_point(2)
      
! check the required orientation of the line i.e. is end point 1 or end point 2 the same
! point as the last point set in the local_cable_segment_list

        if ( same_cell_point(last_cell_point,line_end_cell_point1) ) then

! set up the parameters for a forward loop over the line segments      
          first_cell_segment=1
	  last_cell_segment=n_line_cell_segments
	  cell_segment_step=1
          first_cell_segment_point =1
          second_cell_segment_point=2
	
          write(cable_info_file_unit,*)'First line point: 1'
      
        else if (same_cell_point(last_cell_point,line_end_cell_point2)) then

! set up the parameters for a backward loop over the line segments      
          first_cell_segment=n_line_cell_segments
	  last_cell_segment=1
	  cell_segment_step=-1
          first_cell_segment_point =2
          second_cell_segment_point=1
	
          write(cable_info_file_unit,*)'First line point: 2'
      
        else
! cable route is not continuous so we have an error
          GOTO 9000
        end if
      
!loop over the line cell segments copying into the local_cable_segment_list
        do cell_segment_number=first_cell_segment,last_cell_segment,cell_segment_step
	
          segment_count=segment_count+1
      
          local_cable_segment_list(segment_count)%segment_point(1)=	&
	        problem_lines(line)%cell_segment_list(cell_segment_number)%segment_point(first_cell_segment_point)
		
          local_cable_segment_list(segment_count)%segment_point(2)=	&
	        problem_lines(line)%cell_segment_list(cell_segment_number)%segment_point(second_cell_segment_point)
	
	  local_direction_sign_list(segment_count)=segment_sign
	
	  segment_sign=segment_sign*(-1)
	
          write(cable_info_file_unit,*)'Setting segment',segment_count,	&
	                             ' direction',local_direction_sign_list(segment_count)
          write(cable_info_file_unit,*)'      ',local_cable_segment_list(segment_count)%segment_point(1)%cell,   &
                                face_string(local_cable_segment_list(segment_count)%segment_point(1)%point)
          write(cable_info_file_unit,*)'      ',local_cable_segment_list(segment_count)%segment_point(2)%cell,   &
                                face_string(local_cable_segment_list(segment_count)%segment_point(2)%point)

          last_cell_point=local_cable_segment_list(segment_count)%segment_point(2)
	  
        end do ! next cell segment in the line
      
      end if ! this line has at least one cell segment on it.
      
    end do ! next line_loop (line on the cable route)
    
! STAGE 7: if end point 2 is on a face then create a segment to bring the 
!          current point to the termination face.  

    if (cable_end_cell_point2%point.NE.centre) then
    
      segment_count=segment_count+1
      local_cable_segment_list(segment_count)%segment_point(1)=	&
            local_cable_segment_list(segment_count-1)%segment_point(2)
      
! ensure that the face point is in the same cell as the cell centre point and swap sides
! if required 
	    
      cell_point1=cable_end_cell_point2
      call get_other_side_of_face(cell_point1,cell_point2)
      
      if ( (cell_point1%cell%i.eq.local_cable_segment_list(segment_count)%segment_point(1)%cell%i).AND. &
           (cell_point1%cell%j.eq.local_cable_segment_list(segment_count)%segment_point(1)%cell%j).AND. &
           (cell_point1%cell%k.eq.local_cable_segment_list(segment_count)%segment_point(1)%cell%k) ) then
      
! set termination segment point to cell_point1	 
        local_cable_segment_list(segment_count)%segment_point(2)=cell_point1
      
      else if ( (cell_point2%cell%i.eq.local_cable_segment_list(segment_count)%segment_point(1)%cell%i).AND. &
                (cell_point2%cell%j.eq.local_cable_segment_list(segment_count)%segment_point(1)%cell%j).AND. &
                (cell_point2%cell%k.eq.local_cable_segment_list(segment_count)%segment_point(1)%cell%k) ) then
		
! set termination segment point to cell_point2 
        local_cable_segment_list(segment_count)%segment_point(2)=cell_point2
      
      else ! error - neither point is in the same cell
      
        last_cell_point=local_cable_segment_list(segment_count)%segment_point(1)
        GOTO 9010
	
      end if 
	
      local_direction_sign_list(segment_count)=-1 ! cable direction is from centre to face

	
      write(cable_info_file_unit,*)'Setting segment',segment_count,' direction=',	&
                                   local_direction_sign_list(segment_count)
      write(cable_info_file_unit,*)'      ',local_cable_segment_list(segment_count)%segment_point(1)%cell,   &
                                face_string(local_cable_segment_list(segment_count)%segment_point(1)%point)
      write(cable_info_file_unit,*)'      ',local_cable_segment_list(segment_count)%segment_point(2)%cell,   &
                                face_string(local_cable_segment_list(segment_count)%segment_point(2)%point)
  
    end if
    
    write(cable_info_file_unit,*)'Initial number of segments=',segment_count
    write(cable_info_file_unit,*)'Filter cable segment list'
    
    CALL filter_cable_segment_list(local_cable_segment_list,local_direction_sign_list,	&
                                   local_number_of_cable_segments,segment_count)

! write cable segment list to the cable info file unit

    write(cable_info_file_unit,*)'Final number of segments=',segment_count
    write(cable_info_file_unit,*)'Final cable segment list'

   do i=1,segment_count
	
      write(cable_info_file_unit,*)'Set segment',i,' direction=',	&
                                   local_direction_sign_list(i)
      write(cable_info_file_unit,*)'      ',local_cable_segment_list(i)%segment_point(1)%cell,   &
                                face_string(local_cable_segment_list(i)%segment_point(1)%point)
      write(cable_info_file_unit,*)'      ',local_cable_segment_list(i)%segment_point(2)%cell,   &
                                face_string(local_cable_segment_list(i)%segment_point(2)%point)
   
   end do ! next cable segment

! we now have the complete cell_segment list for the cable so save to the cable_list structure

    cable_list(cable)%number_of_cable_segments=segment_count
    ALLOCATE( cable_list(cable)%cable_segment_list(1:segment_count) )
    cable_list(cable)%cable_segment_list(1:segment_count)=local_cable_segment_list(1:segment_count)
    
! allocate the direction sign list     
    ALLOCATE( cable_list(cable)%direction_sign_list(1:segment_count) )
    cable_list(cable)%direction_sign_list(1:segment_count)=local_direction_sign_list(1:segment_count)

! allocate the bundle segment list     
    ALLOCATE( cable_list(cable)%bundle_segment_list(1:segment_count) )
    cable_list(cable)%bundle_segment_list(1:segment_count)=0
    
    DEALLOCATE( local_cable_segment_list )    
    
    DEALLOCATE( local_direction_sign_list )    
 
  end do ! next cable

  CALL write_line('FINISHED: create_cable_meshes',0,output_to_screen_flag)
    
  RETURN
  
9000 CALL write_line('Error in create_cable_meshes:',0,.TRUE.)
     CALL write_line('Line cell segment list is not continuous',0,.TRUE.)
     write(*,*)'Line=',line
     write(*,*)'Last cell point:',last_cell_point
     write(*,*)'Line end point1:',line_end_cell_point1
     write(*,*)'Line end point2:',line_end_cell_point2
     STOP
  
9010 CALL write_line('Error in create_cable_meshes:',0,.TRUE.)
     CALL write_line('Line cell segment list is not continuous',0,.TRUE.)
     write(*,*)'Line=',line
     write(*,*)'Last cell point:',last_cell_point
     write(*,*)'Line end point1:',cell_point1
     write(*,*)'Line end point2:',cell_point2
     STOP
  
9020 CALL write_line('Error in create_cable_meshes:',0,.TRUE.)
     CALL write_line('Cable junction does not exist',0,.TRUE.)
     write(*,*)'Cable=',cable
     write(*,*)'junction1  =',junction1
     write(*,*)'n_cable_junctions=',n_cable_junctions
     STOP
  
9030 CALL write_line('Error in create_cable_meshes:',0,.TRUE.)
     CALL write_line('Cable junction does not exist',0,.TRUE.)
     write(*,*)'Cable=',cable
     write(*,*)'junction2  =',junction2
     write(*,*)'n_cable_junctions=',n_cable_junctions
     STOP
  
END SUBROUTINE create_cable_meshes

!SUBROUTINE filter_cable_segment_list
!
! NAME
!     SUBROUTINE filter_cable_segment_list
!
! DESCRIPTION
!      Filter out doubled up cell segments at cable terminations if the exist
! for example...
!
! Setting segment          10  direction           1
!                 11          11          11   zmin
!                 11          11          11 centre
! Setting segment          11  direction=          -1
!                 11          11          11 centre
!                 11          11          11   zmin
! 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 28/11/2012 CJS
!
!
SUBROUTINE filter_cable_segment_list(segment_list,direction_sign_list,	&
                                     list_length,n_segments)

USE Geometry_types
USE cell_parameters

IMPLICIT NONE

  integer		:: list_length
  type(cell_segment)	:: segment_list(1:list_length)
  integer		:: direction_sign_list(1:list_length)
  integer		:: n_segments

! local variables

  integer		:: n_segments_in
  integer		:: new_segment,old_segment
  integer		:: i
  
  type(cell_segment)	:: new_segment_list(1:list_length)
  integer		:: new_direction_sign_list(1:list_length)
  
  type(cell_segment)	:: null_segment
  
! function variables  
  logical	:: same_cell_point
  
! START

  n_segments_in=n_segments
  
  if (n_segments_in.le.2) RETURN

  null_segment%segment_point(1)%cell%i=0
  null_segment%segment_point(1)%cell%j=0
  null_segment%segment_point(1)%cell%k=0
  null_segment%segment_point(1)%point =0
  null_segment%segment_point(2)%cell%i=0
  null_segment%segment_point(2)%cell%j=0
  null_segment%segment_point(2)%cell%k=0
  null_segment%segment_point(2)%point =0

! reset the new lists  
  do new_segment=1,list_length
    new_segment_list(new_segment)=null_segment
    new_direction_sign_list(new_segment)=0
  end do

! copy segments into the new_segment list as required

  new_segment=0
  old_segment=0
  
! check end 1 termination
  if ( .NOT.same_cell_point(segment_list(1)%segment_point(1),	&
                            segment_list(2)%segment_point(2)) ) then
! copy the first two segments
 
    old_segment=old_segment+1
    new_segment=new_segment+1
    new_segment_list(new_segment)=segment_list(old_segment)
    new_direction_sign_list(new_segment)=direction_sign_list(old_segment)
    old_segment=old_segment+1
    new_segment=new_segment+1
    new_segment_list(new_segment)=segment_list(old_segment)
    new_direction_sign_list(new_segment)=direction_sign_list(old_segment)

  else
! skip the first two segments

    old_segment=old_segment+2
    
  end if

! we are now at segment 3: loop over segment list copying the segments until we reach the other end
  do i=3,n_segments_in-2
  
    old_segment=old_segment+1
    new_segment=new_segment+1
    new_segment_list(new_segment)=segment_list(old_segment) 
    new_direction_sign_list(new_segment)=direction_sign_list(old_segment)
    
  end do
  
! check end 2 termination
  
  if ( .NOT.same_cell_point(segment_list(n_segments_in-1)%segment_point(1),	&
                            segment_list(n_segments_in)%segment_point(2)) ) then
! copy the last two segments
 
    old_segment=old_segment+1
    new_segment=new_segment+1
    new_segment_list(new_segment)=segment_list(old_segment)
    new_direction_sign_list(new_segment)=direction_sign_list(old_segment)
    old_segment=old_segment+1
    new_segment=new_segment+1
    new_segment_list(new_segment)=segment_list(old_segment)
    new_direction_sign_list(new_segment)=direction_sign_list(old_segment)

  else
! skip the first two segments

    old_segment=old_segment+2
    
  end if
  
! copy the new segment list over the old segment list and return
  n_segments=new_segment
  segment_list(1:list_length)=new_segment_list(1:list_length)
  direction_sign_list(1:list_length)=new_direction_sign_list(1:list_length)
  
  RETURN
  
END SUBROUTINE filter_cable_segment_list

!SUBROUTINE get_closest_mesh_face
!
! NAME
!     SUBROUTINE get_closest_mesh_face
!
! DESCRIPTION
!      Return the closest face to the specified mesh point which has a face set in the mesh
!      This is used to terminate wires to surfaces properly
! 
!     
! COMMENTS
!     re-ordering the loops will improve efficiency enormously...
!
! HISTORY
!
!     started 13/2/2013 CJS
!             17/6/2013 Loop over surfaces with material properties set only - fixes a problem where
!                       cables could snap to output surfaces.
!

  SUBROUTINE  get_closest_mesh_face(point,cable_face)

USE TLM_general
USE Cables
USE Geometry_types
USE Geometry
USE TLM_surface_materials
USE File_information

IMPLICIT NONE
  
  integer	:: point,cable_face

! local variables  
  
  type(ijk)		:: cell
  type(cell_point)	:: cell_face_point
  
  type(xyz)		:: xyz_point
  
  integer		:: surface,face
  integer		:: material_number,i
  integer		:: test_face
  type(cell_point)	:: test_cell_face_point
  
  type(xyz)		:: xyz_face_centre
  real*8		:: distance
  real*8		:: min_dist
  integer		:: min_face
  integer		:: min_surface
  logical		:: face_found

! function types  

  logical		:: same_cell_point
  real*8		:: xyz_distance
  
! START  

! cell containing the point
  cell=problem_points(point)%cell
  cell_face_point%cell=cell

! coordinates of the point  
  xyz_point=problem_points(point)%point
  
  write(*,*)'CALLED:get_closest_mesh_face'
  write(*,*)'point=',point
  write(*,*)'Coordinates',xyz_point%x,xyz_point%y,xyz_point%z
  write(*,*)'Cell ',cell_face_point%cell%i,cell_face_point%cell%j,cell_face_point%cell%k
  
  min_dist=1d30
  min_face=0
  min_surface=0
  face_found=.FALSE.
  
! loop over the faces in the cell  
  do test_face=1,6
  
    cell_face_point%point=test_face
    
! check whether this cell face is in the mesh and has a material allocated to it.

! loop over surface materials
    do material_number=1,n_surface_materials
    
! loop over the geometric surfaces with this material type 
      do i=1,surface_material_list(material_number)%n_surfaces
    
        surface=surface_material_list(material_number)%surface_list(i)
      
!    do surface=1,n_surfaces ! old loop over surfaces - causes problems as output surfaces are also checked...

        do face=1,problem_surfaces(surface)%number_of_faces
      
          test_cell_face_point=problem_surfaces(surface)%face_list(face)
	
	  if (same_cell_point(cell_face_point,test_cell_face_point)) then
	
            face_found=.TRUE.

! get the coordinates of the found face centre
	    CALL get_cell_point_coordinate(test_cell_face_point,xyz_face_centre)
	    	  
! calculate the distance between the found face centre and the termination point
	    distance=xyz_distance(xyz_point,xyz_face_centre)
	  
	    if (distance.lt.min_dist) then
	      min_dist=distance
	      min_face=test_face
	      min_surface=surface
	    end if
	
	  end if
	
        end do ! next face
	
      end do ! next geometric surface with with this material type
      
    end do ! next surface material

! if this face is in the mesh, see if it is the closest to the specified point  
  
  end do ! next face in the cell  
  
  if (face_found) then
  
    cable_face=min_face  
    
    write(cable_info_file_unit,*)'CALLED:get_closest_mesh_face'
    write(cable_info_file_unit,*)'point=',point,' cable_face=',cable_face
    write(cable_info_file_unit,*)'Coordinates',xyz_point%x,xyz_point%y,xyz_point%z
    write(cable_info_file_unit,*)'Cell ',cell_face_point%cell%i,cell_face_point%cell%j,cell_face_point%cell%k
    write(cable_info_file_unit,*)'Minimum distance=',min_dist
    write(cable_info_file_unit,*)'surface number=',min_surface
    write(cable_info_file_unit,*)'face=',min_face

  else
  
    write(*,*)'Error in get_closest_mesh_face'
    write(*,*)'No mesh surface found for cable termination'
    write(*,*)'Point number',point
    write(*,*)'Coordinates',xyz_point%x,xyz_point%y,xyz_point%z
    write(*,*)'Cell ',cell_face_point%cell%i,cell_face_point%cell%j,cell_face_point%cell%k
    STOP
    
  end if
  
  RETURN

  END SUBROUTINE  get_closest_mesh_face
