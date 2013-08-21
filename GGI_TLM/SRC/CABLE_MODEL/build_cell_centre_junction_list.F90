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
!SUBROUTINE build_cell_centre_junction_list
!
! NAME
!     SUBROUTINE build_cell_centre_junction_list
!
! DESCRIPTION
!       construct the data structure for the solver to use
!       in its cell by cell update loop through the mesh
!       The code treats every cell centre as a junction 
!     
! This includes junctions defined in the input file and those formed by the 
! cables traversing the cells along their route.
! The bundle junctions must include the possibility of more than one defined junction and
! also multiple cables within a single cell.
!
! The process is as follows:
!
! STAGE 1: LOOP OVER BUNDLE SEGMENTS, PUTTING THEM INTO THE MESH
! STAGE 2: LOOP OVER THE MESH COUNTING THE SEGMENT CELLS
! STAGE 3: ALLOCATE AND RESET CELL_CENTRE_JUNCTION LIST
! STAGE 4: SET CELL AND BUNDLE SEGMENT DATA IN THE CELL_CENTRE_JUNCTION_LIST
! STAGE 5: SET DEFINED CABLE JUNCTIONS IN THE CELL_CENTRE_JUNCTION_LIST
! STAGE 6: SET THE CABLE_LIST DATA RELATING TO JUNCTIONS
! STAGE 7: ADD THE STRAIGHT THROUGH CABLE CONNECTION NODES TO THE CELL JUNCTIONS
! STAGE 8: ALLOCATE P MATRICES FOR EACH OF THE CELL CENTRE JUNCTIONS
! STAGE 9: LOOP OVER CABLES ADDING TERMINATION JUNCTIONS AND MID-CABLE CONNECTIONS 
!
!     
! COMMENTS
!     This is quite complicated, messy and potentially causes problems. 
!     Ideally we should think of a neater way to do this 
!
! STAGE 7 question:
!    cable_list(cable)%n_conductors=n_conductors   ! Don't know why we do this here - maybe it should be done earlier?
!
! HISTORY
!
!     started 19/09/2012 CJS
!             17/12/2012 CJS	Started to implement junction loads
!     frequency warping included 12/02/2013 CJS
!
!
SUBROUTINE build_cell_centre_junction_list()

USE TLM_general
USE Cables
USE geometry
USE cell_parameters
USE filter_types
USE filter_functions

IMPLICIT NONE

! local variables

  integer,allocatable	:: bundle_segment_number(:,:,:,:)
  integer,allocatable	:: segment_count(:,:,:)
  integer,allocatable	:: cell_to_bundle_junction(:,:,:)
  
  integer		:: n_cable_cells
  
  type(cell_segment)	:: cable_segment
  integer		:: ix1,iy1,iz1,face1
  integer		:: ix2,iy2,iz2,face2
  integer		:: ix,iy,iz,face

  integer		:: segment,segment2
  integer		:: first_segment,last_segment
  integer		:: cable_cell_loop
  integer		:: cable_cell
  
  integer		:: cable_geometry
  
  integer		:: n_internal,n_external_face
  integer		:: n_conductors,conductor
  
  integer		:: junction
  integer		:: point_number

  integer		:: cable
  integer		:: n

  integer		:: bundle_junction
  integer		:: first_internal
  integer		:: first_external
  integer		:: n_external
  integer		:: i
  integer		:: point_in_list
  
  integer		:: impedance_filter
  integer		:: internal_node1,internal_node2
  
  TYPE(Zfilter) :: Zfilter1
  TYPE(Zfilter) :: Zfilter2
  real*8	:: Z_f
  
  integer	:: filter_aorder
  integer	:: filter_border
  
! START

  CALL write_line('CALLED: build_cell_centre_junction_list',0,output_to_screen_flag)
  
  ALLOCATE( bundle_segment_number(1:nx,1:ny,1:nz,1:6) )
  ALLOCATE( segment_count(1:nx,1:ny,1:nz) )
  ALLOCATE( cell_to_bundle_junction(1:nx,1:ny,1:nz) )

  bundle_segment_number(1:nx,1:ny,1:nz,1:6)=0
  segment_count(1:nx,1:ny,1:nz)=0
  cell_to_bundle_junction(1:nx,1:ny,1:nz)=0
  
! STAGE 1: LOOP OVER BUNDLE SEGMENTS, PUTTING THEM INTO THE MESH
  do segment=1,n_bundle_segments

! get segment coordinates  
    cable_segment=bundle_segment_list(segment)%cable_segment
    ix1  =cable_segment%segment_point(1)%cell%i
    iy1  =cable_segment%segment_point(1)%cell%j
    iz1  =cable_segment%segment_point(1)%cell%k
    face1=cable_segment%segment_point(1)%point
    ix2  =cable_segment%segment_point(2)%cell%i
    iy2  =cable_segment%segment_point(2)%cell%j
    iz2  =cable_segment%segment_point(2)%cell%k
    face2=cable_segment%segment_point(2)%point

! check that the points are in the same cell
    if ( (ix1.ne.ix2).OR.(iy1.ne.iy2).OR.(iz1.ne.iz2) ) then
      GOTO 9000
    end if
    
    if (face1.ne.centre) then
      face=face1
    else if (face2.ne.centre) then
      face=face2
    else
      GOTO 9010    
    end if
    
    segment_count(ix1,iy1,iz1)=segment_count(ix1,iy1,iz1)+1
    
    if (bundle_segment_number(ix1,iy1,iz1,face).EQ.0) then
      bundle_segment_number(ix1,iy1,iz1,face)=segment
    else
      GOTO 9020
    end if
    
  end do ! next bundle segment

! STAGE 2: LOOP OVER THE MESH COUNTING THE SEGMENT CELLS
  n_cable_cells=0
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx

        if (segment_count(ix,iy,iz).NE.0) then
	  n_cable_cells=n_cable_cells+1
	end if

      end do 
    end do 
  end do 
  
  CALL write_line_integer('Number of cable cells=',n_cable_cells,0,output_to_screen_flag)
  
! STAGE 3: ALLOCATE AND RESET CELL_CENTRE_JUNCTION LIST

  n_cell_centre_junctions=n_cable_cells
  ALLOCATE( cell_centre_junction_list(1:n_cable_cells) )

! reset cell centre junction list
  do cable_cell_loop=1,n_cable_cells
    cell_centre_junction_list(cable_cell_loop)%n_segments=6
    ALLOCATE( cell_centre_junction_list(cable_cell_loop)%segment_list(1:6) )
    ALLOCATE( cell_centre_junction_list(cable_cell_loop)%n_external_conductors(1:7) )
    ALLOCATE( cell_centre_junction_list(cable_cell_loop)%external_conductor_count(1:7) )
    ALLOCATE( cell_centre_junction_list(cable_cell_loop)%P_matrix_list(1:7) )
    cell_centre_junction_list(cable_cell_loop)%n_internal_connection_nodes=0
    cell_centre_junction_list(cable_cell_loop)%internal_connection_node_count=0
    cell_centre_junction_list(cable_cell_loop)%n_external_conductors(1:7)=0
    cell_centre_junction_list(cable_cell_loop)%external_conductor_count(1:7)=0
    cell_centre_junction_list(cable_cell_loop)%segment_list(1:6)=0
    cell_centre_junction_list(cable_cell_loop)%number_of_cable_junctions=0
    cell_centre_junction_list(cable_cell_loop)%n_internal_impedance_filters=0
  end do
  
! STAGE 4: SET CELL AND BUNDLE SEGMENT DATA IN THE CELL_CENTRE_JUNCTION_LIST

  n_cable_cells=0
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx

        if (segment_count(ix,iy,iz).NE.0) then
	
	  n_cable_cells=n_cable_cells+1
	  
	  cell_centre_junction_list(n_cable_cells)%cell_point%cell%i=ix
	  cell_centre_junction_list(n_cable_cells)%cell_point%cell%j=iy
	  cell_centre_junction_list(n_cable_cells)%cell_point%cell%k=iz
	  cell_centre_junction_list(n_cable_cells)%cell_point%point=centre
	  
	  cell_to_bundle_junction(ix,iy,iz)=n_cable_cells
	  
	  do face=1,6
	    segment=bundle_segment_number(ix,iy,iz,face)
	    cell_centre_junction_list(n_cable_cells)%segment_list(face)=segment
	    if (segment.ne.0) then
	      cell_centre_junction_list(n_cable_cells)%n_external_conductors(face)=bundle_segment_list(segment)%n_conductors
	    else
	      cell_centre_junction_list(n_cable_cells)%n_external_conductors(face)=0
	    end if
	  end do
	  
	end if

      end do 
    end do 
  end do 

! STAGE 5: SET DEFINED CABLE JUNCTIONS IN THE CELL_CENTRE_JUNCTION_LIST
! first count the number of cable junctions in the cell centre

  do junction=1,n_cable_junctions
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then

! Already set in create_cable_meshes...    
!      point_number=cable_junction_list(junction)%point_number
!      cable_junction_list(junction)%cell_point%cell=problem_points(point_number)%cell

      ix=cable_junction_list(junction)%cell_point%cell%i
      iy=cable_junction_list(junction)%cell_point%cell%j
      iz=cable_junction_list(junction)%cell_point%cell%k
    
      cable_cell=cell_to_bundle_junction(ix,iy,iz)
    
      if (cable_cell.ne.0) then
    
! set cable cell (cell centre junction) data 
        cell_centre_junction_list(cable_cell)%number_of_cable_junctions=	&
             cell_centre_junction_list(cable_cell)%number_of_cable_junctions+1
	   
        cell_centre_junction_list(cable_cell)%internal_connection_node_count=	&
            cell_centre_junction_list(cable_cell)%internal_connection_node_count+	&
	    cable_junction_list(junction)%n_internal_connection_nodes
	  
        cable_junction_list(junction)%bundle_junction_number=cable_cell

! ********** count the number of internal impedances associated with this cable_junction **********
        cell_centre_junction_list(cable_cell)%n_internal_impedance_filters=	&
            cell_centre_junction_list(cable_cell)%n_internal_impedance_filters+	&
	    cable_junction_list(junction)%number_of_internal_impedances

        cell_centre_junction_list(cable_cell)%n_external_conductors(7)=	&
            cell_centre_junction_list(cable_cell)%n_external_conductors(7)+	&
	    cable_junction_list(junction)%number_of_internal_impedances
      
      else
! there should be a cell centre junction defined, if not we have an error...
        GOTO 9030
      end if  
    
    end if !  junction_type.EQ.junction_type_cell
      
  end do ! next cable junction

! allocate the cell_centre_junction_list
  do cable_cell=1,n_cable_cells
  
    n=cell_centre_junction_list(cable_cell)%number_of_cable_junctions
    if (n.gt.0) then
      ALLOCATE( cell_centre_junction_list(cable_cell)%cable_junction_list(1:n) )
    end if
    
! reset number_of_cable_junctions so it can be used as a counter to fill the cable_junction_list 
    
    cell_centre_junction_list(cable_cell)%number_of_cable_junctions=0
    
  end do

! fill cell_centre_junction_list(cable_cell)%cable_junction_list
 
  do junction=1,n_cable_junctions
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
  
      cable_cell=cable_junction_list(junction)%bundle_junction_number
    
      cell_centre_junction_list(cable_cell)%number_of_cable_junctions=	&
           cell_centre_junction_list(cable_cell)%number_of_cable_junctions+1
	     
      n=cell_centre_junction_list(cable_cell)%number_of_cable_junctions
      cell_centre_junction_list(cable_cell)%cable_junction_list(n)=junction
    
      cable_junction_list(junction)%first_internal_connection_node_in_bundle_junction=	&
          cell_centre_junction_list(cable_cell)%n_internal_connection_nodes+1
    
      cell_centre_junction_list(cable_cell)%n_internal_connection_nodes=	&
                  cell_centre_junction_list(cable_cell)%n_internal_connection_nodes+	&
	  	cable_junction_list(junction)%n_internal_connection_nodes
    
    end if !  junction_type.EQ.junction_type_cell
    
  end do ! next cable junction

! STAGE 6: SET THE CABLE_LIST DATA RELATING TO JUNCTIONS

  do cable=1,n_cables

! junction 1  
    junction=cable_list(cable)%junction_1
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
    
      cable_list(cable)%bundle_junction_1=cable_junction_list(junction)%bundle_junction_number
    
      cable_list(cable)%first_internal_connection_node_in_bundle_junction_1=	&
            cable_junction_list(junction)%first_internal_connection_node_in_bundle_junction
      cable_list(cable)%first_external_conductor_in_bundle_junction_1=0

    end if !  junction_type.EQ.junction_type_cell
    
! junction 2
    junction=cable_list(cable)%junction_2
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
    
      cable_list(cable)%bundle_junction_2=cable_junction_list(junction)%bundle_junction_number
      cable_list(cable)%first_internal_connection_node_in_bundle_junction_2=	&
            cable_junction_list(junction)%first_internal_connection_node_in_bundle_junction
      cable_list(cable)%first_external_conductor_in_bundle_junction_2=0

    end if !  junction_type.EQ.junction_type_cell
  
  end do ! next cable
  
! STAGE 7: ADD THE STRAIGHT THROUGH CABLE CONNECTION NODES TO THE CELL JUNCTIONS
  do cable=1,n_cables
  
    cable_geometry=cable_list(cable)%cable_geometry_number
    n_conductors=cable_geometry_list(cable_geometry)%n_conductors
    cable_list(cable)%n_conductors=n_conductors   ! Don't know why we do this here - maybe it should be done earlier?

! get the first and last segment numbers for the cell centre segment loop    
    junction=cable_list(cable)%junction_1
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
      first_segment=2
    else
      first_segment=1
    end if !  junction_type.EQ.junction_type_cell
    
    junction=cable_list(cable)%junction_2
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
      last_segment=cable_list(cable)%number_of_cable_segments-1
    else
      last_segment=cable_list(cable)%number_of_cable_segments
    end if !  junction_type.EQ.junction_type_cell
    
! loop over internal segments (note step 2 because each cell has two segments associated with it)

    do segment=first_segment,last_segment,2
    
      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)
      
      cable_cell=cell_to_bundle_junction(ix,iy,iz)

! each cable contributes a number of internal connection nodes equal to cable_list(cable)%n_conductors 

      cell_centre_junction_list(cable_cell)%n_internal_connection_nodes=	&
           cell_centre_junction_list(cable_cell)%n_internal_connection_nodes+n_conductors		  
    
    end do ! next internal segment on cable
  
  end do ! next cable

! We now have the number of internal connection nodes and external conductors for each of the cable cells
! so we can allocate the P matrices 

! STAGE 8: ALLOCATE P MATRICES FOR EACH OF THE CELL CENTRE JUNCTIONS
  
  do cable_cell_loop=1,n_cell_centre_junctions
    
    n_internal=cell_centre_junction_list(cable_cell_loop)%n_internal_connection_nodes
    
    do face=1,7 ! include internal impedances on face 7.
    
      n_external_face=cell_centre_junction_list(cable_cell_loop)%n_external_conductors(face)
      if (n_external_face.ne.0) then
      
        ALLOCATE( cell_centre_junction_list(cable_cell_loop)%P_matrix_list(face)%P(1:n_internal,1:n_external_face) )

! reset P matrix for this cell face
        cell_centre_junction_list(cable_cell_loop)%P_matrix_list(face)%P(1:n_internal,1:n_external_face)=0	
	
      end if 
      
    end do ! next face
    
  end do ! next cable_cell  

! STAGE 9: LOOP OVER CABLES ADDING END POINT JUNCTIONS AND MID-CABLE CONNECTIONS 
! to the cell_centre_junction_list P matrices

! ********** First allocate internal impedance filter data ********** 

  do cable_cell=1,n_cable_cells
  
    n=cell_centre_junction_list(cable_cell)%n_internal_impedance_filters
    if (n.gt.0) then
      ALLOCATE( cell_centre_junction_list(cable_cell)%Sfilter(1:n) )
      ALLOCATE( cell_centre_junction_list(cable_cell)%Zfilter(1:n) )
      ALLOCATE( cell_centre_junction_list(cable_cell)%Z_f(1:n) )
    end if
    
  end do

! note that the following loop mimics that which builds the bundle segments. This
! is necessary in order to get the numbering correct. 
  
  do cable=1,n_cables
  
    n_conductors=cable_list(cable)%n_conductors
    n_external=n_conductors
    
! End 1 junction
    segment=1
    junction=cable_list(cable)%junction_1
  
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
    
      bundle_junction=cable_list(cable)%bundle_junction_1
      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)
      cable_cell=cell_to_bundle_junction(ix,iy,iz)
    
! check consistency here:
      if (cable_cell.ne.bundle_junction) then
        write(*,*)'Inconsistent bundle junction number'
        write(*,*)'Cable=',cable
        write(*,*)'End 1 junction',junction
        write(*,*)'End 1 bundle junction',bundle_junction
        write(*,*)'Cable_cell=',cable_cell
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
      
      if(point_in_list.eq.0) then
        write(*,*)'Error in create_bundle_junctions'
        write(*,*)'cable end not found in junction cable list'
        write(*,*)'cable=',cable,' end=1'
        write(*,*)'junction=',junction
        write(*,*)'number_of_cables=',cable_junction_list(junction)%number_of_cables
	write(*,*)
        do i=1,cable_junction_list(junction)%number_of_cables
	  write(*,*)i,cable_junction_list(junction)%cable_list(i),cable_junction_list(junction)%cable_end_list(i)
	end do
        STOP
      end if
    
      first_external=cell_centre_junction_list(cable_cell)%external_conductor_count(face)+1

! copy P matrix from cable_junction_list into cell_centre_junction_list   

      cell_centre_junction_list(cable_cell)%P_matrix_list(face)%	&
           P(first_internal:first_internal+n_internal-1,first_external:first_external+n_external-1)= &
           cable_junction_list(junction)%Pmatrix(point_in_list)%P(1:n_internal,1:n_external) 

! add the conductors to the external conductor count      
      cell_centre_junction_list(cable_cell)%external_conductor_count(face)=	&
           cell_centre_junction_list(cable_cell)%external_conductor_count(face)+n_external

! ************* Add the internal impedance data for this junction to the cell_centre_junction_list  ***********

      do impedance_filter=1,cable_junction_list(junction)%number_of_internal_impedances

! Set P matrix
        internal_node1=cable_junction_list(junction)%node_1(impedance_filter)
        internal_node2=cable_junction_list(junction)%node_2(impedance_filter)
	
	if (internal_node1.ne.0) then
          cell_centre_junction_list(cable_cell)%P_matrix_list(7)%	&
             P(first_internal+internal_node1-1,first_external+impedance_filter-1)=1
        end if
	if (internal_node2.ne.0) then
          cell_centre_junction_list(cable_cell)%P_matrix_list(7)%	&
             P(first_internal+internal_node2-1,first_external+impedance_filter-1)=-1
        end if
	
! Copy Sfilter
	cell_centre_junction_list(cable_cell)%Sfilter(first_external+impedance_filter-1)=	&
	        cable_junction_list(junction)%Sfilter(impedance_filter)
	
! Calculate Zfilter
!        Zfilter1=s_to_z(cable_junction_list(junction)%Sfilter(impedance_filter),dt) 
        Zfilter1=s_to_z_warp(cable_junction_list(junction)%Sfilter(impedance_filter),dt,bicubic_warp_flag,frequency_scale) 
        CALL Z_fast_slow_docomposition( Zfilter1 ,Z_f  , Zfilter2 )
	
	cell_centre_junction_list(cable_cell)%Zfilter(first_external+impedance_filter-1)=Zfilter2
	cell_centre_junction_list(cable_cell)%Z_f(first_external+impedance_filter-1)=Z_f

      end do ! next impedance_filter
	     
    end if ! junction1.EQ.cell_centre_junction

! get the first and last segment numbers for the cell centre segment loop    
    junction=cable_list(cable)%junction_1
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
      first_segment=2
    else
      first_segment=1
    end if !  junction_type.EQ.junction_type_cell
    
    junction=cable_list(cable)%junction_2
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
      last_segment=cable_list(cable)%number_of_cable_segments-1
    else
      last_segment=cable_list(cable)%number_of_cable_segments
    end if !  junction_type.EQ.junction_type_cell
      
! loop over internal segments (note step 2 because each cell has two segments associated with it)
     
    do segment=first_segment,last_segment,2

      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k
      cable_cell=cell_to_bundle_junction(ix,iy,iz)

! each cable contributes a number of internal connection nodes equal to cable_list(cable)%n_conductors 

      first_internal=cell_centre_junction_list(cable_cell)%internal_connection_node_count+1

! add the conductors to the internal conductor count      
      cell_centre_junction_list(cable_cell)%internal_connection_node_count=	&
           cell_centre_junction_list(cable_cell)%internal_connection_node_count+n_conductors
     
      if (cell_centre_junction_list(cable_cell)%internal_connection_node_count.gt.	&
          cell_centre_junction_list(cable_cell)%n_internal_connection_nodes) then
	  
	write(*,*)'Internal connection node counting error'
	write(*,*)'count=',cell_centre_junction_list(cable_cell)%internal_connection_node_count
	write(*,*)'total=',cell_centre_junction_list(cable_cell)%n_internal_connection_nodes
	STOP
	  
      end if

! do first segment of the pair through the cell    
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)
      first_external=cell_centre_junction_list(cable_cell)%external_conductor_count(face)+1

! Set P matrix elements for segment 1    
      do i=1,n_conductors
        cell_centre_junction_list(cable_cell)%P_matrix_list(face)%P(first_internal+i-1,first_external+i-1)=1
      end do

! add the conductors to the external conductor count      
      cell_centre_junction_list(cable_cell)%external_conductor_count(face)=	&
           cell_centre_junction_list(cable_cell)%external_conductor_count(face)+n_conductors

! do second segment of the pair through the cell    
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment+1),face)
      first_external=cell_centre_junction_list(cable_cell)%external_conductor_count(face)+1

! Set P matrix elements for segment 1    
      do i=1,n_conductors
        cell_centre_junction_list(cable_cell)%P_matrix_list(face)%P(first_internal+i-1,first_external+i-1)=1
      end do

! add the conductors to the external conductor count      
      cell_centre_junction_list(cable_cell)%external_conductor_count(face)=	&
           cell_centre_junction_list(cable_cell)%external_conductor_count(face)+n_conductors
    
    end do ! next internal segment on cable
    
    
! End 2 junction
    segment=cable_list(cable)%number_of_cable_segments
    junction=cable_list(cable)%junction_2
    
    if (cable_junction_list(junction)%junction_type.EQ.junction_type_cell) then
    
      bundle_junction=cable_list(cable)%bundle_junction_2
      ix=cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%i
      iy=cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%j
      iz=cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%k
      CALL get_cell_segment_face(cable_list(cable)%cable_segment_list(segment),face)
      cable_cell=cell_to_bundle_junction(ix,iy,iz)
    
! check consistency here:
      if (cable_cell.ne.bundle_junction) then
        write(*,*)'Inconsistent bundle junction number'
	write(*,*)'This may occur if the junction position and '
	write(*,*)'the end of the cable are not at the same point in space'
        write(*,*)'Cable=',cable
        write(*,*)'End 2 junction',junction
        write(*,*)'End 2 bundle junction',bundle_junction
        write(*,*)'Cable_cell=',cable_cell
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
      
      if(point_in_list.eq.0) then
        write(*,*)'Error in create_bundle_junctions'
        write(*,*)'cable end not found in junction cable list'
        write(*,*)'cable=',cable,' end=2'
        write(*,*)'junction=',junction
        write(*,*)'number_of_cables=',cable_junction_list(junction)%number_of_cables
	write(*,*)
        do i=1,cable_junction_list(junction)%number_of_cables
	  write(*,*)i,cable_junction_list(junction)%cable_list(i),cable_junction_list(junction)%cable_end_list(i)
	end do
        STOP
      end if
     
      first_external=cell_centre_junction_list(cable_cell)%external_conductor_count(face)+1

! copy P matrix from cable_junction_list into cell_centre_junction_list   

      cell_centre_junction_list(cable_cell)%P_matrix_list(face)%	&
           P(first_internal:first_internal+n_internal-1,first_external:first_external+n_external-1)= &
           cable_junction_list(junction)%Pmatrix(point_in_list)%P(1:n_internal,1:n_external)

! add the conductors to the external conductor count      
      cell_centre_junction_list(cable_cell)%external_conductor_count(face)=	&
           cell_centre_junction_list(cable_cell)%external_conductor_count(face)+n_external

! ************* Add the internal impedance data for this junction to the cell_centre_junction_list  ***********

      do impedance_filter=1,cable_junction_list(junction)%number_of_internal_impedances

! Set P matrix
        internal_node1=cable_junction_list(junction)%node_1(impedance_filter)
        internal_node2=cable_junction_list(junction)%node_2(impedance_filter)
	
	if (internal_node1.ne.0) then
          cell_centre_junction_list(cable_cell)%P_matrix_list(7)%	&
             P(first_internal+internal_node1-1,first_external+impedance_filter-1)=1
        end if
	if (internal_node2.ne.0) then
          cell_centre_junction_list(cable_cell)%P_matrix_list(7)%	&
             P(first_internal+internal_node2-1,first_external+impedance_filter-1)=-1
        end if
	
! Copy Sfilter
	cell_centre_junction_list(cable_cell)%Sfilter(first_external+impedance_filter-1)=	&
	        cable_junction_list(junction)%Sfilter(impedance_filter)
	
! Calculate Zfilter
! Calculate Z domain susceptibility admittance filter coefficients by bilinear transformation
!        Zfilter1=s_to_z(cable_junction_list(junction)%Sfilter(impedance_filter),dt) 
        Zfilter1=s_to_z_warp(cable_junction_list(junction)%Sfilter(impedance_filter),dt,bicubic_warp_flag,frequency_scale) 
        CALL Z_fast_slow_docomposition( Zfilter1 ,Z_f  , Zfilter2 )
	
	cell_centre_junction_list(cable_cell)%Zfilter(first_external+impedance_filter-1)=Zfilter2
	cell_centre_junction_list(cable_cell)%Z_f(first_external+impedance_filter-1)=Z_f

      end do ! next impedance_filter
  
    end if ! junction.EQ.junction_type_cell
  
  end do ! next cable
  
  DEALLOCATE( bundle_segment_number )
  DEALLOCATE( segment_count )
  DEALLOCATE( cell_to_bundle_junction )

  CALL write_line('FINISHED: build_cell_centre_junction_list',0,output_to_screen_flag)
    
  RETURN
  
9000 CALL write_line('ERROR in build_cell_centre_junction_list',0,.TRUE.)
     CALL write_line('Cable segment points are not in the same cell',0,.TRUE.)
     CALL write_line_integer('Segment',segment,0,.TRUE.)
     STOP
     
9010 CALL write_line('ERROR in build_cell_centre_junction_list',0,.TRUE.)
     CALL write_line('No face point in cell segment',0,.TRUE.)
     CALL write_line_integer('Segment',segment,0,.TRUE.)
     STOP
     
9020 CALL write_line('ERROR in build_cell_centre_junction_list',0,.TRUE.)
     CALL write_line('cell segment already set',0,.TRUE.)
     CALL write_line_integer('Segment',segment,0,.TRUE.)
     CALL write_line_integer('ix=',ix1,0,.TRUE.)
     CALL write_line_integer('iy=',iy1,0,.TRUE.)
     CALL write_line_integer('iz=',iz1,0,.TRUE.)
     CALL write_line_integer('face=',face,0,.TRUE.)
     STOP
     
9030 CALL write_line('ERROR in build_cell_centre_junction_list',0,.TRUE.)
     CALL write_line('No bundle junction defined for the cable junction',0,.TRUE.)
     CALL write_line_integer('Cable junction number',junction,0,.TRUE.)
     STOP
     
     STOP
  
END SUBROUTINE build_cell_centre_junction_list
