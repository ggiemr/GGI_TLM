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
!SUBROUTINE build_face_junction_list
!
! NAME
!     SUBROUTINE build_face_junction_list
!
! DESCRIPTION
!       construct the data structure for the solver to use
!       in its face by face update connect loop through the mesh
!       The face_junction_list includes cable terminations to surfaces
!
!       The process is similar to that used to build the cell_centre_junction_list
!
!
! The process is as follows:
!
! STAGE 1: LOOP OVER BUNDLE SEGMENTS, PUTTING THEM INTO THE MESH
! STAGE 2: LOOP OVER THE MESH COUNTING THE CABLE FACES
! STAGE 3: ALLOCATE AND RESET FACE_JUNCTION LIST
! STAGE 4: SET CELL AND BUNDLE SEGMENT DATA IN THE FACE_JUNCTION_LIST
! STAGE 5: SET DEFINED CABLE JUNCTIONS IN THE FACE_JUNCTION_LIST
! STAGE 6: SET THE CABLE_LIST DATA RELATING TO JUNCTIONS
! STAGE 7: ADD THE STRAIGHT THROUGH CABLE CONNECTION NODES TO THE FACE JUNCTIONS
! STAGE 8: ALLOCATE P MATRICES FOR EACH OF THE FACE JUNCTIONS
! STAGE 9: LOOP OVER CABLES ADDING TERMINATIONS, BOUNDARY CONDITIONS AND MID-CABLE CONNECTIONS 
!
! COMMENTS
!     
!
! HISTORY
!
!     started 21/09/2012 CJS
!     frequency warping included 12/02/2013 CJS
!
!
SUBROUTINE build_face_junction_list()

USE TLM_general
USE Cables
USE geometry
USE cell_parameters
USE filter_types
USE filter_functions

IMPLICIT NONE

! local variables

  integer,allocatable	:: bundle_segment_number(:,:,:,:)
  integer,allocatable	:: face_to_bundle_junction(:,:,:,:)
  
  integer		:: n_cable_faces
  
  integer		:: cable_geometry
  
  type(cell_segment)	:: cable_segment
  integer		:: ix1,iy1,iz1,face1
  integer		:: ix2,iy2,iz2,face2
  integer		:: ix,iy,iz,face
  
  integer		:: jface

  integer		:: segment
  integer		:: cable_face_loop
  integer		:: cable_face
  
  integer		:: cable
  integer		:: point_number
  
  integer		:: n_internal,n_external_face
  integer		:: n_conductors,conductor
  
  integer		:: junction
  integer		:: bundle_junction
  
  integer		:: first_internal
  integer		:: first_external
  integer		:: first_segment
  integer		:: last_segment
  integer		:: n_external
  integer		:: i,n
  integer		:: point_in_list
  
  integer		:: impedance_filter
  integer		:: internal_node1,internal_node2
  
  TYPE(Zfilter) :: Zfilter1
  TYPE(Zfilter) :: Zfilter2
  real*8	:: Z_f
  
  integer	:: filter_aorder
  integer	:: filter_border
 
! START

  CALL write_line('CALLED: build_face_junction_list',0,output_to_screen_flag)
  
  ALLOCATE( bundle_segment_number(1:nx,1:ny,1:nz,1:6) )
  ALLOCATE( face_to_bundle_junction(1:nx,1:ny,1:nz,1:6) )

  bundle_segment_number(1:nx,1:ny,1:nz,1:6)=0
  face_to_bundle_junction(1:nx,1:ny,1:nz,1:6)=0
  
! STAGE 1: LOOP OVER BUNDLE SEGMENTS, PUTTING THEM INTO THE MESH
  
  do segment=1,n_bundle_segments

! get segment coordinates  
    cable_segment=bundle_segment_list(segment)%cable_segment
    ix  =cable_segment%segment_point(1)%cell%i
    iy  =cable_segment%segment_point(1)%cell%j
    iz  =cable_segment%segment_point(1)%cell%k
    CALL get_cell_segment_face(cable_segment%segment_point,face)      
    
    if (bundle_segment_number(ix,iy,iz,face).EQ.0) then
      bundle_segment_number(ix,iy,iz,face)=segment
    else
      GOTO 9020
    end if
    
  end do ! next bundle segment

! STAGE 2: LOOP OVER THE MESH COUNTING THE CABLE FACES
  
  n_cable_faces=0

! loop over the xmin faces of the mesh looking for cables intersecting faces
  do iz=1,nz
    do iy=1,ny
      do ix=2,nx

        if ( (bundle_segment_number(ix  ,iy,iz,face_xmin).NE.0).OR.	&
	     (bundle_segment_number(ix-1,iy,iz,face_xmax).NE.0) ) then
! a cable intersects this face
	  n_cable_faces=n_cable_faces+1
	end if

      end do 
    end do 
  end do 

! loop over the ymin faces of the mesh looking for cables intersecting faces
  do iz=1,nz
    do iy=2,ny
      do ix=1,nx

        if ( (bundle_segment_number(ix,iy  ,iz,face_ymin).NE.0).OR.	&
	     (bundle_segment_number(ix,iy-1,iz,face_ymax).NE.0) ) then
! a cable intersects this face
	  n_cable_faces=n_cable_faces+1
	end if

      end do 
    end do 
  end do 

! loop over the zmin faces of the mesh looking for cables intersecting faces
  do iz=2,nz
    do iy=1,ny
      do ix=1,nx

        if ( (bundle_segment_number(ix,iy,iz  ,face_zmin).NE.0).OR.	&
	     (bundle_segment_number(ix,iy,iz-1,face_zmax).NE.0) ) then
! a cable intersects this face
	  n_cable_faces=n_cable_faces+1
	end if

      end do 
    end do 
  end do 
  
  CALL write_line_integer('Number of cable faces=',n_cable_faces,0,output_to_screen_flag)
  
! STAGE 3: ALLOCATE AND RESET FACE_JUNCTION LIST

  n_face_junctions=n_cable_faces
  ALLOCATE( face_junction_list(1:n_face_junctions) )

! reset cable_face bundle segments

  do cable_face_loop=1,n_face_junctions
    face_junction_list(cable_face_loop)%n_segments=2
    ALLOCATE( face_junction_list(cable_face_loop)%segment_list(1:2) )
    ALLOCATE( face_junction_list(cable_face_loop)%n_external_conductors(1:3) )
    ALLOCATE( face_junction_list(cable_face_loop)%external_conductor_count(1:3) )
    ALLOCATE( face_junction_list(cable_face_loop)%P_matrix_list(1:3) )
    face_junction_list(cable_face_loop)%n_internal_connection_nodes=0
    face_junction_list(cable_face_loop)%internal_connection_node_count=0
    face_junction_list(cable_face_loop)%number_of_cable_junctions=0
    face_junction_list(cable_face_loop)%segment_list(1:2)=0
    face_junction_list(cable_face_loop)%n_external_conductors(1:3)=0    
    face_junction_list(cable_face_loop)%external_conductor_count(1:3)=0    
    face_junction_list(cable_face_loop)%n_internal_impedance_filters=0
  end do
  
! STAGE 4: SET CELL AND BUNDLE SEGMENT DATA IN THE FACE_JUNCTION_LIST

  n_cable_faces=0

! loop over the xmin faces of the mesh looking for cables intersecting faces
  do iz=1,nz
    do iy=1,ny
      do ix=2,nx

        if ( (bundle_segment_number(ix  ,iy,iz,face_xmin).NE.0).OR.	&
	     (bundle_segment_number(ix-1,iy,iz,face_xmax).NE.0) ) then

! a cable intersects this face
	  n_cable_faces=n_cable_faces+1
	  
	  face_junction_list(n_cable_faces)%cell_point%cell%i=ix
	  face_junction_list(n_cable_faces)%cell_point%cell%j=iy
	  face_junction_list(n_cable_faces)%cell_point%cell%k=iz
	  face_junction_list(n_cable_faces)%cell_point%point=face_xmin
	  
	  face_to_bundle_junction(ix  ,iy,iz,face_xmin)=n_cable_faces	  	  	  	  
	  face_to_bundle_junction(ix-1,iy,iz,face_xmax)=n_cable_faces	  	  	  	  
	  	  
	  segment=bundle_segment_number(ix-1,iy,iz,face_xmax)
	  face_junction_list(n_cable_faces)%segment_list(1)=segment
	  if (segment.ne.0) then
	    face_junction_list(n_cable_faces)%n_external_conductors(1)=bundle_segment_list(segment)%n_conductors
	  else
	    face_junction_list(n_cable_faces)%n_external_conductors(1)=0
	  end if
	  	  
	  segment=bundle_segment_number(ix  ,iy,iz,face_xmin)	  
	  face_junction_list(n_cable_faces)%segment_list(2)=segment
	  if (segment.ne.0) then
	    face_junction_list(n_cable_faces)%n_external_conductors(2)=bundle_segment_list(segment)%n_conductors
	  else
	    face_junction_list(n_cable_faces)%n_external_conductors(2)=0
	  end if
	  
	end if

      end do 
    end do 
  end do 

! loop over the ymin faces of the mesh looking for cables intersecting faces
  do iz=1,nz
    do iy=2,ny
      do ix=1,nx

        if ( (bundle_segment_number(ix,iy  ,iz,face_ymin).NE.0).OR.	&
	     (bundle_segment_number(ix,iy-1,iz,face_ymax).NE.0) ) then
! a cable intersects this face

	  n_cable_faces=n_cable_faces+1
	  
	  face_junction_list(n_cable_faces)%cell_point%cell%i=ix
	  face_junction_list(n_cable_faces)%cell_point%cell%j=iy
	  face_junction_list(n_cable_faces)%cell_point%cell%k=iz
	  face_junction_list(n_cable_faces)%cell_point%point=face_ymin
	  
	  face_to_bundle_junction(ix,iy  ,iz,face_ymin)=n_cable_faces	  	  	  	  
	  face_to_bundle_junction(ix,iy-1,iz,face_ymax)=n_cable_faces	  	  	  	  
	  
	  segment=bundle_segment_number(ix,iy-1,iz,face_ymax)	  
	  face_junction_list(n_cable_faces)%segment_list(1)=segment
	  if (segment.ne.0) then
	    face_junction_list(n_cable_faces)%n_external_conductors(1)=bundle_segment_list(segment)%n_conductors
	  else
	    face_junction_list(n_cable_faces)%n_external_conductors(1)=0
	  end if
	  
	  segment=bundle_segment_number(ix,iy  ,iz,face_ymin)	  
	  face_junction_list(n_cable_faces)%segment_list(2)=segment
	  if (segment.ne.0) then
	    face_junction_list(n_cable_faces)%n_external_conductors(2)=bundle_segment_list(segment)%n_conductors
	  else
	    face_junction_list(n_cable_faces)%n_external_conductors(2)=0
	  end if

	end if

      end do 
    end do 
  end do 

! loop over the zmin faces of the mesh looking for cables intersecting faces
  do iz=2,nz
    do iy=1,ny
      do ix=1,nx

        if ( (bundle_segment_number(ix,iy,iz  ,face_zmin).NE.0).OR.	&
	     (bundle_segment_number(ix,iy,iz-1,face_zmax).NE.0) ) then
! a cable intersects this face

	  n_cable_faces=n_cable_faces+1
	  
	  face_junction_list(n_cable_faces)%cell_point%cell%i=ix
	  face_junction_list(n_cable_faces)%cell_point%cell%j=iy
	  face_junction_list(n_cable_faces)%cell_point%cell%k=iz
	  face_junction_list(n_cable_faces)%cell_point%point=face_zmin
	  
	  face_to_bundle_junction(ix,iy,iz  ,face_zmin)=n_cable_faces	  	  
	  face_to_bundle_junction(ix,iy,iz-1,face_zmax)=n_cable_faces	  	  
	  
	  segment=bundle_segment_number(ix,iy,iz-1,face_zmax)	  
	  face_junction_list(n_cable_faces)%segment_list(1)=segment
	  if (segment.ne.0) then
	    face_junction_list(n_cable_faces)%n_external_conductors(1)=bundle_segment_list(segment)%n_conductors
	  else
	    face_junction_list(n_cable_faces)%n_external_conductors(1)=0
	  end if
	  
	  segment=bundle_segment_number(ix,iy,iz  ,face_zmin)	  
	  face_junction_list(n_cable_faces)%segment_list(2)=segment
	  if (segment.ne.0) then
	    face_junction_list(n_cable_faces)%n_external_conductors(2)=bundle_segment_list(segment)%n_conductors
	  else
	    face_junction_list(n_cable_faces)%n_external_conductors(2)=0
	  end if

	end if

      end do 
    end do 
  end do 
  
! STAGE 5: SET DEFINED CABLE JUNCTIONS IN THE FACE_JUNCTION_LIST
! first count the number of cable junctions in the cell centre

  do junction=1,n_cable_junctions
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then

! Already set in create_cable_meshes...  
!      point_number=cable_junction_list(junction)%point_number
!      cable_junction_list(junction)%cell_point=problem_points(point_number)%face
      
      ix=cable_junction_list(junction)%cell_point%cell%i
      iy=cable_junction_list(junction)%cell_point%cell%j
      iz=cable_junction_list(junction)%cell_point%cell%k
      face=cable_junction_list(junction)%cell_point%point
      
! note: the junction face may not be a material face here so we need to check and create a warning...
    
      cable_face=face_to_bundle_junction(ix,iy,iz,face)
    
      if (cable_face.ne.0) then
    
! set cable face (face_junction) data 
        face_junction_list(cable_face)%number_of_cable_junctions=	&
             face_junction_list(cable_face)%number_of_cable_junctions+1
	   
        face_junction_list(cable_face)%internal_connection_node_count=	&
            face_junction_list(cable_face)%internal_connection_node_count+	&
	    cable_junction_list(junction)%n_internal_connection_nodes
	  
        cable_junction_list(junction)%bundle_junction_number=cable_face

! ********** count the number of internal impedances associated with this cable_junction **********
        face_junction_list(cable_face)%n_internal_impedance_filters=	&
            face_junction_list(cable_face)%n_internal_impedance_filters+	&
	    cable_junction_list(junction)%number_of_internal_impedances

        face_junction_list(cable_face)%n_external_conductors(3)=	&
            face_junction_list(cable_face)%n_external_conductors(3)+	&
	    cable_junction_list(junction)%number_of_internal_impedances
      
      else
! there should be a face junction defined, if not we have an error...
        GOTO 9030
      end if  
    
    end if !  junction_type.EQ.junction_type_face
      
  end do ! next cable junction

! allocate the face_junction_list
  do cable_face=1,n_face_junctions
  
    n=face_junction_list(cable_face)%number_of_cable_junctions
    if (n.gt.0) then
      ALLOCATE( face_junction_list(cable_face)%cable_junction_list(1:n) )
    end if
    
! reset number_of_cable_junctions so it can be used as a counter to fill the cable_junction_list 
    
    face_junction_list(cable_face)%number_of_cable_junctions=0
    
  end do

! fill face_junction_list(cable_face)%cable_junction_list
 
  do junction=1,n_cable_junctions
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
  
      cable_face=cable_junction_list(junction)%bundle_junction_number
    
      face_junction_list(cable_face)%number_of_cable_junctions=	&
           face_junction_list(cable_face)%number_of_cable_junctions+1
	     
      n=face_junction_list(cable_face)%number_of_cable_junctions
      face_junction_list(cable_face)%cable_junction_list(n)=junction
    
      cable_junction_list(junction)%first_internal_connection_node_in_bundle_junction=	&
          face_junction_list(cable_face)%n_internal_connection_nodes+1
    
      face_junction_list(cable_face)%n_internal_connection_nodes=	&
                  face_junction_list(cable_face)%n_internal_connection_nodes+	&
	  	cable_junction_list(junction)%n_internal_connection_nodes
    
    end if !  junction_type.EQ.junction_type_face
    
  end do ! next cable junction
  
! STAGE 6: SET THE CABLE_LIST DATA RELATING TO JUNCTIONS

  do cable=1,n_cables

! junction 1  
    junction=cable_list(cable)%junction_1
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
    
      cable_list(cable)%bundle_junction_1=cable_junction_list(junction)%bundle_junction_number
    
      cable_list(cable)%first_internal_connection_node_in_bundle_junction_1=	&
            cable_junction_list(junction)%first_internal_connection_node_in_bundle_junction
      cable_list(cable)%first_external_conductor_in_bundle_junction_1=0

    end if !  junction_type.EQ.junction_type_face
    
! junction 2
    junction=cable_list(cable)%junction_2
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
    
      cable_list(cable)%bundle_junction_2=cable_junction_list(junction)%bundle_junction_number
      cable_list(cable)%first_internal_connection_node_in_bundle_junction_2=	&
            cable_junction_list(junction)%first_internal_connection_node_in_bundle_junction
      cable_list(cable)%first_external_conductor_in_bundle_junction_2=0

    end if !  junction_type.EQ.junction_type_face
  
  end do ! next cable

! STAGE 7: ADD THE STRAIGHT THROUGH CABLE CONNECTION NODES TO THE FACE JUNCTIONS
  
  do cable=1,n_cables
  
    cable_geometry=cable_list(cable)%cable_geometry_number
    n_conductors=cable_geometry_list(cable_geometry)%n_conductors

! get the first and last segment numbers for the face segment loop    
    junction=cable_list(cable)%junction_1
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
      first_segment=2
    else! junction1 is at a cell centre
      first_segment=1
    end if !  junction_type.EQ.junction_type_face
    
    junction=cable_list(cable)%junction_2
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
      last_segment=cable_list(cable)%number_of_cable_segments-2
    else  ! junction2 is at a cell centre
      last_segment=cable_list(cable)%number_of_cable_segments-1
    end if !  junction_type.EQ.junction_type_face
    
! loop over internal segments (note step 2 because each face has two segments associated with it)

    do segment=first_segment,last_segment,2
    
      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)      
      
      cable_face=face_to_bundle_junction(ix,iy,iz,face)

! each cable contributes a number of internal connection nodes equal to cable_list(cable)%n_conductors 

      face_junction_list(cable_face)%n_internal_connection_nodes=	&
           face_junction_list(cable_face)%n_internal_connection_nodes+n_conductors		  
    
    end do ! next internal segment on cable
  
  end do ! next cable

! We now have the number of internal connection nodes and external conductors for each of the cable faces
! so we can allocate the P matrices 

! STAGE 8: ALLOCATE P MATRICES and BOUNDARY CONDITION VECTOR FOR EACH OF THE FACE JUNCTIONS

  do cable_face_loop=1,n_face_junctions
    
    n_internal=face_junction_list(cable_face_loop)%n_internal_connection_nodes
    
    ALLOCATE( face_junction_list(cable_face_loop)%BC(1:n_internal) )
    
! reset BC vector for this cell face

    face_junction_list(cable_face_loop)%BC(1:n_internal)=0	
    
    do face=1,3
    
      n_external_face=face_junction_list(cable_face_loop)%n_external_conductors(face)

      if (n_external_face.ne.0) then
      		  
        ALLOCATE( face_junction_list(cable_face_loop)%P_matrix_list(face)%P(1:n_internal,1:n_external_face) )

! reset P matrix for this cell face

        face_junction_list(cable_face_loop)%P_matrix_list(face)%P(1:n_internal,1:n_external_face)=0	
	
      end if 
      
    end do ! next face
    
  end do ! next cable_face  

! STAGE 9: LOOP OVER CABLES ADDING TERMINATIONS AND MID-CABLE CONNECTIONS 
! TO THE FACE_JUNCTION_LIST P MATRICES
! ALSO SET BOUNDARY CONDITIONS

! note that the following loop mimics that which builds the bundle segments. This
! is necessary in order to get the numbering correct. 

! ********** First allocate internal impedance filter data ********** 
  do cable_face=1,n_cable_faces
  
    n=face_junction_list(cable_face)%n_internal_impedance_filters
    if (n.gt.0) then
      ALLOCATE( face_junction_list(cable_face)%Sfilter(1:n) )
      ALLOCATE( face_junction_list(cable_face)%Zfilter(1:n) )
      ALLOCATE( face_junction_list(cable_face)%Z_f(1:n) )
    end if
    
  end do
  
  do cable=1,n_cables
  
    n_conductors=cable_list(cable)%n_conductors
    n_external=n_conductors
    
    write(*,*)'cable',cable,' n_external=',n_external
    
! End 1 junction
    
    segment=1
    junction=cable_list(cable)%junction_1
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
      
      bundle_junction=cable_list(cable)%bundle_junction_1
      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)      
      cable_face=face_to_bundle_junction(ix,iy,iz,face)
      
      if ( (face.eq.face_xmin).OR.(face.eq.face_ymin).OR.(face.eq.face_zmin) )then
        jface=2
      else
        jface=1
      end if
    
! check consistency here:
      if (cable_face.ne.bundle_junction) then
        write(*,*)'Inconsistent bundle junction number'
        write(*,*)'Cable=',cable
        write(*,*)'End 1 junction',junction
        write(*,*)'End 1 bundle junction',bundle_junction
        write(*,*)'cable_face=',cable_face
        STOP
      end if

      first_internal=cable_list(cable)%first_internal_connection_node_in_bundle_junction_1
      n_internal=cable_junction_list(junction)%n_internal_connection_nodes

! find the current cable and end in the cable_junction(junction)%cable list

      point_in_list=0
      do i=1,cable_junction_list(junction)%number_of_cables
        if ( (cable_junction_list(junction)%cable_list(i).eq.cable).AND.	&
             (cable_junction_list(junction)%cable_end_list(i).eq.1) ) then
	   
	  if(point_in_list.ne.0) then
	    write(*,*)'Error in create_bundle_junctions'
	    write(*,*)'cable end found more than once in junction cable list'
	    write(*,*)'cable=',cable
	    write(*,*)'junction=',junction
	    STOP
	  end if
	
          point_in_list=i
	
        end if
      end do ! next cable in junction cable list
    
      first_external=face_junction_list(cable_face)%external_conductor_count(jface)+1

! copy P matrix from cable_junction_list into face_junction_list   
      
      face_junction_list(cable_face)%P_matrix_list(jface)%	&
           P(first_internal:first_internal+n_internal-1,first_external:first_external+n_external-1)= &
           cable_junction_list(junction)%Pmatrix(point_in_list)%P(1:n_internal,1:n_external) 

! copy BC vector from cable_junction_list into face_junction_list   
      face_junction_list(cable_face)%BC(first_internal:first_internal+n_internal-1)= &
           cable_junction_list(junction)%BC(1:n_internal) 

! add the conductors to the external conductor count      
      face_junction_list(cable_face)%external_conductor_count(jface)=	&
           face_junction_list(cable_face)%external_conductor_count(jface)+n_external

! ************* Add the internal impedance data for this junction to the cell_centre_junction_list  ***********

      do impedance_filter=1,cable_junction_list(junction)%number_of_internal_impedances

! Set P matrix
        internal_node1=cable_junction_list(junction)%node_1(impedance_filter)
        internal_node2=cable_junction_list(junction)%node_2(impedance_filter)
	if (internal_node1.ne.0) then
          face_junction_list(cable_face)%P_matrix_list(3)%	&
             P(first_internal+internal_node1-1,first_external+impedance_filter-1)=1
        end if
	if (internal_node2.ne.0) then
          face_junction_list(cable_face)%P_matrix_list(3)%	&
             P(first_internal+internal_node2-1,first_external+impedance_filter-1)=-1
        end if
	
! Copy Sfilter
	face_junction_list(cable_face)%Sfilter(first_external+impedance_filter-1)=	&
	        cable_junction_list(junction)%Sfilter(impedance_filter)
	
! Calculate Zfilter
! Calculate Z domain susceptibility admittance filter coefficients by bilinear transformation
!        Zfilter1=s_to_z(cable_junction_list(junction)%Sfilter(impedance_filter),dt) 
        Zfilter1=s_to_z_warp(cable_junction_list(junction)%Sfilter(impedance_filter),dt,bicubic_warp_flag,frequency_scale) 
        CALL Z_fast_slow_docomposition( Zfilter1 ,Z_f  , Zfilter2 )
	
	face_junction_list(cable_face)%Zfilter(first_external+impedance_filter-1)=Zfilter2
	face_junction_list(cable_face)%Z_f(first_external+impedance_filter-1)=Z_f
	
      end do ! next impedance_filter
	     
    end if ! junction1.EQ.face_junction

! get the first and last segment numbers for the face segment loop    
    junction=cable_list(cable)%junction_1
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
      first_segment=2
    else! junction1 is at a cell centre
      first_segment=1
    end if !  junction_type.EQ.junction_type_face
    
    junction=cable_list(cable)%junction_2
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
      last_segment=cable_list(cable)%number_of_cable_segments-2
    else  ! junction2 is at a cell centre
      last_segment=cable_list(cable)%number_of_cable_segments-1
    end if !  junction_type.EQ.junction_type_face
      
! loop over internal segments (note step 2 because each face has two segments associated with it)
     
    do segment=first_segment,last_segment,2

      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)      
      cable_face=face_to_bundle_junction(ix,iy,iz,face)

! each cable contributes a number of internal connection nodes equal to cable_list(cable)%n_conductors 

      first_internal=face_junction_list(cable_face)%internal_connection_node_count+1

! add the conductors to the internal conductor count      
      face_junction_list(cable_face)%internal_connection_node_count=	&
           face_junction_list(cable_face)%internal_connection_node_count+n_conductors
     
      if (face_junction_list(cable_face)%internal_connection_node_count.gt.	&
          face_junction_list(cable_face)%n_internal_connection_nodes) then
	  
	write(*,*)'Internal connection node counting error'
	write(*,*)'count=',face_junction_list(cable_face)%internal_connection_node_count
	write(*,*)'total=',face_junction_list(cable_face)%n_internal_connection_nodes
	STOP
	  
      end if

! do first segment of the pair through the face    
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)      
      if ( (face.eq.face_xmin).OR.(face.eq.face_ymin).OR.(face.eq.face_zmin) )then
        jface=2
      else
        jface=1
      end if
      first_external=face_junction_list(cable_face)%external_conductor_count(jface)+1

! Set P matrix elements for segment 1    
      do i=1,n_conductors
        face_junction_list(cable_face)%P_matrix_list(jface)%P(first_internal+i-1,first_external+i-1)=1
      end do

! add the conductors to the external conductor count      
      face_junction_list(cable_face)%external_conductor_count(jface)=	&
           face_junction_list(cable_face)%external_conductor_count(jface)+n_conductors

! do second segment of the pair through the face   
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment+1),face)
      if ( (face.eq.face_xmin).OR.(face.eq.face_ymin).OR.(face.eq.face_zmin) )then
        jface=2
      else
        jface=1
      end if
      first_external=face_junction_list(cable_face)%external_conductor_count(jface)+1

! Set P matrix elements for segment 2  
      do i=1,n_conductors
        face_junction_list(cable_face)%P_matrix_list(jface)%P(first_internal+i-1,first_external+i-1)=1
      end do

! add the conductors to the external conductor count      
      face_junction_list(cable_face)%external_conductor_count(jface)=	&
           face_junction_list(cable_face)%external_conductor_count(jface)+n_conductors
    
    end do ! next internal segment on cable
       
! End 2 junction
    
    segment=cable_list(cable)%number_of_cable_segments
    junction=cable_list(cable)%junction_2
    
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_face) then
    
      bundle_junction=cable_list(cable)%bundle_junction_2
      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)      
      cable_face=face_to_bundle_junction(ix,iy,iz,face)
      
      if ( (face.eq.face_xmin).OR.(face.eq.face_ymin).OR.(face.eq.face_zmin) )then
        jface=2
      else
        jface=1
      end if
    
! check consistency here:
      if (cable_face.ne.bundle_junction) then
        write(*,*)'Inconsistent bundle junction number'
        write(*,*)'Cable=',cable
        write(*,*)'End 2 junction',junction
        write(*,*)'End 2 bundle junction',bundle_junction
        write(*,*)'cable_face=',cable_face
        STOP
      end if

      first_internal=cable_list(cable)%first_internal_connection_node_in_bundle_junction_2
      n_internal=cable_junction_list(junction)%n_internal_connection_nodes

! find the current cable and end in the cable_junction(junction)%cable list

      point_in_list=0
      do i=1,cable_junction_list(junction)%number_of_cables
        if ( (cable_junction_list(junction)%cable_list(i).eq.cable).AND.	&
           (cable_junction_list(junction)%cable_end_list(i).eq.2) ) then
	   
	  if(point_in_list.ne.0) then
	    write(*,*)'Error in create_bundle_junctions'
	    write(*,*)'cable end found more than once in junction cable list'
	    write(*,*)'cable=',cable
	    write(*,*)'junction=',junction
	    STOP
	  end if
	
          point_in_list=i
	
        end if
      end do ! next cable in junction cable list
    
      first_external=face_junction_list(cable_face)%external_conductor_count(jface)+1

! copy P matrix from cable_junction_list into face_junction_list   

      face_junction_list(cable_face)%P_matrix_list(jface)%	&
           P(first_internal:first_internal+n_internal-1,first_external:first_external+n_external-1)= &
           cable_junction_list(junction)%Pmatrix(point_in_list)%P(1:n_internal,1:n_external)

! copy BC vector from cable_junction_list into face_junction_list   
      face_junction_list(cable_face)%BC(first_internal:first_internal+n_internal-1)= &
           cable_junction_list(junction)%BC(1:n_internal) 

! add the conductors to the external conductor count      
      face_junction_list(cable_face)%external_conductor_count(jface)=	&
           face_junction_list(cable_face)%external_conductor_count(jface)+n_external

! ************* Add the internal impedance data for this junction to the cell_centre_junction_list  ***********

      do impedance_filter=1,cable_junction_list(junction)%number_of_internal_impedances

! Set P matrix
        internal_node1=cable_junction_list(junction)%node_1(impedance_filter)
        internal_node2=cable_junction_list(junction)%node_2(impedance_filter)
	if (internal_node1.ne.0) then
          face_junction_list(cable_face)%P_matrix_list(3)%	&
             P(first_internal+internal_node1-1,first_external+impedance_filter-1)=1
        end if
	if (internal_node2.ne.0) then
          face_junction_list(cable_face)%P_matrix_list(3)%	&
             P(first_internal+internal_node2-1,first_external+impedance_filter-1)=-1
        end if
	
! Copy Sfilter
	face_junction_list(cable_face)%Sfilter(first_external+impedance_filter-1)=	&
	        cable_junction_list(junction)%Sfilter(impedance_filter)
	
! Calculate Zfilter
! Calculate Z domain susceptibility admittance filter coefficients by bilinear transformation
!        Zfilter1=s_to_z(cable_junction_list(junction)%Sfilter(impedance_filter),dt) 
        Zfilter1=s_to_z_warp(cable_junction_list(junction)%Sfilter(impedance_filter),dt,bicubic_warp_flag,frequency_scale) 
        CALL Z_fast_slow_docomposition( Zfilter1 ,Z_f  , Zfilter2 )
	
	face_junction_list(cable_face)%Zfilter(first_external+impedance_filter-1)=Zfilter2
	face_junction_list(cable_face)%Z_f(first_external+impedance_filter-1)=Z_f
	
      end do ! next impedance_filter
  
    end if ! junction.EQ.junction_type_face
  
  end do ! next cable
  
  DEALLOCATE( bundle_segment_number )
  DEALLOCATE( face_to_bundle_junction )

  CALL write_line('FINISHED: build_face_junction_list',0,output_to_screen_flag)

  RETURN
  
9000 CALL write_line('ERROR in build_face_junction_list',0,.TRUE.)
     CALL write_line('Cable segment points are not in the same cell',0,.TRUE.)
     CALL write_line_integer('Segment',segment,0,.TRUE.)
     STOP
     
9010 CALL write_line('ERROR in build_face_junction_list',0,.TRUE.)
     CALL write_line('No face point in cell segment',0,.TRUE.)
     CALL write_line_integer('Segment',segment,0,.TRUE.)
     STOP
     
9020 CALL write_line('ERROR in build_face_junction_list',0,.TRUE.)
     CALL write_line('cell segment already set',0,.TRUE.)
     CALL write_line_integer('Segment',segment,0,.TRUE.)
     CALL write_line_integer('ix=',ix1,0,.TRUE.)
     CALL write_line_integer('iy=',iy1,0,.TRUE.)
     CALL write_line_integer('iz=',iz1,0,.TRUE.)
     CALL write_line_integer('face=',face,0,.TRUE.)
     STOP 
     
9030 CALL write_line('ERROR in build_face_junction_list',0,.TRUE.)
     CALL write_line('No bundle junction defined for the cable junction',0,.TRUE.)
     CALL write_line_integer('Cable junction number',junction,0,.TRUE.)
     write(*,*)'Cell coordinates:',ix,iy,iz,' face: ',face_string(face)
     STOP
     
END SUBROUTINE build_face_junction_list
