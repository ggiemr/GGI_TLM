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
! SUBROUTINE set_cable_bundles_in_mesh
! SUBROUTINE initialise_cable_bundles
!
! NAME
!     set_cable_bundles_in_mesh
!
! DESCRIPTION
!    1.  loop over the cell_centre_junction_list and flag in local_cell_cable(i,j,k) 
!    2.  loop over the face_junction_list and flag in local_surface_cable(i,j,k) 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 21/09/2012 CJS
!
!
SUBROUTINE set_cable_bundles_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE Cables
USE file_information

IMPLICIT NONE

! local variables

  integer	:: cx,cy,cz,face
  integer 	:: cell_junction,face_junction
  type(cell_point)	:: cable_face
  
  integer	:: loop
  
! START

  
!  CALL write_line('CALLED: set_cable_bundles_in_mesh',0,output_to_screen_flag)

! Stage 1. Set cell junctions

  do cell_junction=1,n_cell_centre_junctions
  
    cx=cell_centre_junction_list(cell_junction)%cell_point%cell%i
    cy=cell_centre_junction_list(cell_junction)%cell_point%cell%j
    cz=cell_centre_junction_list(cell_junction)%cell_point%cell%k
  
    if (rank.eq.cell_rank(cz)) then
    
      local_cell_cable(cx,cy,cz)=cell_junction
    
      write(*,*)'Setting cell centre junction',cell_junction
      write(*,*)'Coordinates:',cx,cy,cz
    end if
    
  end do ! next cell_centre_junction


! Stage 2. Set face junctions

  do face_junction=1,n_face_junctions
     	     	
    CALL get_min_face(face_junction_list(face_junction)%cell_point,cable_face)
	
    cx=cable_face%cell%i
    cy=cable_face%cell%j
    cz=cable_face%cell%k
    face=cable_face%point
    
    if (rank.eq.cell_face_rank(cz,face)) then
    
      local_surface_cable(cx  ,cy  ,cz  ,face)=face_junction
    
!      write(*,*)'Setting cell face cable point',face_junction
!      write(*,*)'Coordinates:',cx,cy,cz,' face:',face

! Put a check here for proper connection to surfaces
      do loop=1,face_junction_list(face_junction)%n_internal_connection_nodes
        if(face_junction_list(face_junction)%BC(loop).EQ.-1) then
! we have a termination to a surface here - check that a surface exists...

          if (local_surface_material(cx  ,cy  ,cz  ,face).EQ.0) then
	  
	    if (rank.eq.0) then
	      write(warning_file_unit,*)'WARNING: No surface found for cable-surface connection:'
	      write(warning_file_unit,*)'cx=',cx
	      write(warning_file_unit,*)'cy=',cy
	      write(warning_file_unit,*)'cz=',cz
	      write(warning_file_unit,*)'face=',face
	    end if
	    
	    number_of_warnings=number_of_warnings+1
	  
	  end if
	  
        end if
      end do

    end if
    
  end do ! next face_junction  
  
!  CALL write_line('FINISHED: set_cable_bundles_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_cable_bundles_in_mesh
!
! NAME
!     initialise_cable_bundles
!
! DESCRIPTION
!     
! Set cable link and stub matrices for the TLM cable update     
!
! COMMENTS
!     
! Stub calculation: - we need to decide what we need to save on the segment and what is just
! a local array required for the calculation...
!
!
! HISTORY
!
!     started 21/09/2012 CJS
!             17/12/2012 Allocate filter data for internal impedance filters at junctions
!
!
SUBROUTINE initialise_cable_bundles

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE Cables
USE file_information
USE filter_types
USE filter_functions

IMPLICIT NONE

! local variables

  integer	:: segment
  integer	:: n_conductors
  
  real*8,allocatable	:: L(:,:)
  real*8,allocatable	:: C(:,:)
  real*8,allocatable	:: Ci(:,:)
  real*8,allocatable	:: Rs(:,:)
  
  real*8,allocatable	:: Ztl(:,:)
  real*8,allocatable	:: Ytl(:,:)
  
  real*8,allocatable	:: Ztlc(:,:)
  real*8,allocatable	:: Ytlc(:,:)
  
  real*8,allocatable	:: ZLstub(:,:)
  real*8,allocatable	:: ZCstub(:,:)
  real*8,allocatable	:: YCstub(:,:)
  
  real*8,allocatable	:: Zf(:,:)
  real*8,allocatable	:: Yf(:,:)
  
  real*8		:: dt2

  logical	:: write_to_cable_info_file_flag
  
  integer	:: row,col
  integer	:: filter_data
  integer	:: filter
  integer	:: filter_aorder,filter_border
  
  integer	:: cable_cell,cable_face
  integer	:: impedance_filter
  integer	:: n

! START

  CALL write_line('CALLED: initialise_cable_bundles',0,output_to_screen_flag)

  dt2=dt/2d0

  do segment=1,n_bundle_segments
  
    n_conductors=bundle_segment_list(segment)%n_conductors

! Allocate local variables for calculation

    ALLOCATE( L(1:n_conductors,1:n_conductors) )
    ALLOCATE( C(1:n_conductors,1:n_conductors) )
    ALLOCATE( Ci(1:n_conductors,1:n_conductors) )
    ALLOCATE( Rs(1:n_conductors,1:n_conductors) )

    ALLOCATE( Ztl(1:n_conductors,1:n_conductors) )
    ALLOCATE( Ytl(1:n_conductors,1:n_conductors) )
    
    ALLOCATE( ZLstub(1:n_conductors,1:n_conductors) )
    
    ALLOCATE( Zf(1:n_conductors,1:n_conductors) )
    ALLOCATE( Yf(1:n_conductors,1:n_conductors) )

! copy L,C,Rs to local matrices    
    L(1:n_conductors,1:n_conductors)=bundle_segment_list(segment)%L(1:n_conductors,1:n_conductors)
    C(1:n_conductors,1:n_conductors)=bundle_segment_list(segment)%C(1:n_conductors,1:n_conductors)
    Rs(1:n_conductors,1:n_conductors)=bundle_segment_list(segment)%R(1:n_conductors,1:n_conductors)

! Stub calculation 

    CALL calculate_TLM_cable_matrices(L,C,Ztl,Zlstub,n_conductors)
    
    write_to_cable_info_file_flag=.FALSE.
    CALL TLM_Cable_checks(Ztl,ZLstub,n_conductors,write_to_cable_info_file_flag)

    call dsvd_invert(Ztl,n_conductors,n_conductors,Ytl,n_conductors)
        
    Zf(:,:)=Ztl(:,:)+ZLstub(:,:)+Rs(:,:) 

! include fast impedance from filter function into Zf, also count the number of filter responses required
    bundle_segment_list(segment)%n_filter_data=0
    do row=1,n_conductors
      do col=1,n_conductors
      
        filter=bundle_segment_list(segment)%filter_number(row,col)
	if (filter.ne.0) then ! there is a filter on this matrix element
	
	  Zf(row,col)=Zf(row,col)+bundle_segment_list(segment)%Z_f(filter) 
	  bundle_segment_list(segment)%n_filter_data=bundle_segment_list(segment)%n_filter_data+1
	  
	end if
      end do
    end do
			
    call dsvd_invert(Zf,n_conductors,n_conductors,Yf,n_conductors) 

! Allocate bundle_segment matrices

    ALLOCATE( bundle_segment_list(segment)%Zlink(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%Ylink(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%ZLstub(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%Yf(1:n_conductors,1:n_conductors ) )
    
! copy TLM data from local matrices to the bundle_segment structure
    bundle_segment_list(segment)%Zlink(:,:) =Ztl(:,:)
    bundle_segment_list(segment)%Ylink(:,:) =Ytl(:,:)
    bundle_segment_list(segment)%ZLstub(:,:)=ZLstub(:,:)
    bundle_segment_list(segment)%Yf(:,:)=Yf(:,:)
     
! Allocate and reset filter data
    if (bundle_segment_list(segment)%n_filter_data.NE.0) then

      ALLOCATE( bundle_segment_list(segment)%Zfilter_data(1:bundle_segment_list(segment)%n_filter_data) )

! loop over filters allocating memory for the particular filter order
      filter_data=0
      do row=1,n_conductors
        do col=1,n_conductors
      
          filter=bundle_segment_list(segment)%filter_number(row,col)
	  if (filter.ne.0) then ! there is a filter on this matrix element
	
	    filter_data=filter_data+1
	    filter_aorder=bundle_segment_list(segment)%Zfilter(filter)%a%order
	    filter_border=bundle_segment_list(segment)%Zfilter(filter)%b%order
	    
	    bundle_segment_list(segment)%Zfilter_data(filter_data)=	&
	        allocate_Zfilter_response(filter_aorder,filter_border)
	  
	  end if
        end do
      end do
 
    end if  ! n_filter_data.NE.0

! Allocate and reset voltage vectors and current vector
 
    ALLOCATE( bundle_segment_list(segment)%Vlink(1:n_conductors) )
    ALLOCATE( bundle_segment_list(segment)%VLstub(1:n_conductors) )
    ALLOCATE( bundle_segment_list(segment)%Vsource(1:n_conductors) )
    ALLOCATE( bundle_segment_list(segment)%Iw_centre(1:n_conductors) )
    ALLOCATE( bundle_segment_list(segment)%Iw_face(1:n_conductors) )
    
    bundle_segment_list(segment)%Vlink(1:n_conductors)=0d0
    bundle_segment_list(segment)%VLstub(1:n_conductors)=0d0
    bundle_segment_list(segment)%Vsource(1:n_conductors)=0d0
    bundle_segment_list(segment)%Iw_centre(1:n_conductors)=0d0
    bundle_segment_list(segment)%Iw_face(1:n_conductors)=0d0

! Deallocate local variables for calculation

    DEALLOCATE( L )
    DEALLOCATE( C )
    DEALLOCATE( Ci )
    DEALLOCATE( Rs )

    DEALLOCATE( Ztl )
    DEALLOCATE( Ytl )
    
    DEALLOCATE( ZLstub )
    
    DEALLOCATE( Zf )
    DEALLOCATE( Yf )

  end do ! next segment
  
! Allocate and set data relating to cell centre junction internal impedances
  do cable_cell=1,n_cell_centre_junctions
  
    n=cell_centre_junction_list(cable_cell)%n_internal_impedance_filters
    if (n.gt.0) then
      ALLOCATE( cell_centre_junction_list(cable_cell)%Yf(1:n,1:n) )
      cell_centre_junction_list(cable_cell)%Yf(1:n,1:n)=0d0
      ALLOCATE( cell_centre_junction_list(cable_cell)%Zfilter_data(1:n) )
    end if

    do impedance_filter=1,cell_centre_junction_list(cable_cell)%n_internal_impedance_filters
	
! Set Yf(:,:)
      cell_centre_junction_list(cable_cell)%Yf(impedance_filter,impedance_filter)=	&
           1d0/cell_centre_junction_list(cable_cell)%Z_f(impedance_filter)
		
! Allocate Zfilter_data
      filter_aorder=cell_centre_junction_list(cable_cell)%Zfilter(impedance_filter)%a%order
      filter_border=cell_centre_junction_list(cable_cell)%Zfilter(impedance_filter)%b%order
      
      cell_centre_junction_list(cable_cell)%Zfilter_data(impedance_filter)=     &
          allocate_Zfilter_response(filter_aorder,filter_border)

    end do ! next impedance_filter
    
  end do ! next cell centre junction
  
! Allocate and set data relating to face centre junction internal impedances
  do cable_face=1,n_face_junctions
  
    n=face_junction_list(cable_face)%n_internal_impedance_filters
    if (n.gt.0) then
      ALLOCATE( face_junction_list(cable_face)%Yf(1:n,1:n) )
      face_junction_list(cable_face)%Yf(1:n,1:n)=0d0
      ALLOCATE( face_junction_list(cable_face)%Zfilter_data(1:n) )
    end if

    do impedance_filter=1,face_junction_list(cable_face)%n_internal_impedance_filters
	
! Set Yf(:,:)
      face_junction_list(cable_face)%Yf(impedance_filter,impedance_filter)=	&
           1d0/face_junction_list(cable_face)%Z_f(impedance_filter)
		
! Allocate Zfilter_data
      filter_aorder=face_junction_list(cable_face)%Zfilter(impedance_filter)%a%order
      filter_border=face_junction_list(cable_face)%Zfilter(impedance_filter)%b%order

      face_junction_list(cable_face)%Zfilter_data(impedance_filter)=     &
          allocate_Zfilter_response(filter_aorder,filter_border)

    end do ! next impedance_filter
    
  end do ! next cell centre junction 
  
  CALL write_line('FINISHED: initialise_cable_bundles',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_cable_bundles
