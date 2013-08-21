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
!SUBROUTINE create_bundle_LCRG_matrices
!
! NAME
!     SUBROUTINE create_bundle_LCRG_matrices
!
! DESCRIPTION
!     calculate the L and C matrices for each bundle segment
!     This subroutine returns the L and C per segment NOT per-unit length
!
!     
! COMMENTS
!     
! In the structure for LC matrix calculation, each cable contributes 
! a single conductor to the external problem at the moment i.e. the shield
! or the single conductor defined in the cable type
!
! HISTORY
!
!     started 19/09/2012 CJS
!     calculation on a bundle segment geometry basis 8/11/2012 CJS
!     include shielded cables 19/11/2012 CJS
!     start to include ribbon cable 24/04/2013 CJS
!
!
SUBROUTINE create_bundle_LCRG_matrices()

USE TLM_general
USE Cables
USE pul_data
USE constants
USE File_information
USE filter_types
USE filter_functions

IMPLICIT NONE

! local variables

  integer	:: bundle_segment
  integer	:: bundle_segment_geometry
  
  integer	:: cable_loop
  integer	:: number_of_cables
  integer	:: cable
  integer	:: cable_geometry
  integer	:: n_conductors
  integer	:: n_shielded_conductors
  integer	:: n_external_conductors
  integer	:: conductor_count,external_conductor_count,total_conductor_count
  
  integer	:: filter_count
  integer	:: n_filters
  integer	:: filter
  
  integer	:: cable_row,cable_col
  integer	:: bundle_row,bundle_col
  integer	:: row,col
  
  integer		:: n_LC_conductors
  integer		:: total_n_bundle_conductors
  
  real*8,allocatable	:: cable_radius(:)
  real*8,allocatable	:: cable_xc(:)
  real*8,allocatable	:: cable_yc(:)
  
  real*8,allocatable	:: conductor_radius(:)
  real*8,allocatable	:: conductor_xc(:)
  real*8,allocatable	:: conductor_yc(:)
  real*8,allocatable	:: dielectric_radius(:)
  real*8,allocatable	:: dielectric_permittivity(:)
  real*8,allocatable	:: local_L(:,:)
  real*8,allocatable	:: local_C(:,:)
  real*8,allocatable	:: Zf(:,:)
  
  integer		:: first_conductor
  integer		:: last_conductor
  integer		:: direction_sign
  integer		:: cable_segment
  integer,allocatable	:: segment_conductor_count(:)
  
  integer,allocatable	:: local_LC_conductor_to_segment_conductor(:)

  real*8	:: segment_length
  
  real*8	:: cable_centre_x,cable_centre_y
  
  real*8	:: reference_radius,reference_radius_L,reference_radius_C

  logical	:: write_to_cable_info_file_flag
  
  integer	:: i
  
! START

  CALL write_line('CALLED: create_bundle_LCRG_matrices',0,output_to_screen_flag)
  
  write(cable_info_file_unit,*)
  write(cable_info_file_unit,*)'create_bundle_LCRG_matrices'
  write(cable_info_file_unit,*)
  
! calculate L,C,R matrices on a bundle geometry basis
  do bundle_segment_geometry=1,n_bundle_segment_geometries
         
    conductor_count=0   
    filter_count=0   
    do cable_loop=1,bundle_segment_geometry_list(bundle_segment_geometry)%n_cables
      
      cable=bundle_segment_geometry_list(bundle_segment_geometry)%cable_list(cable_loop)
      cable_geometry=cable_list(cable)%cable_geometry_number
      n_conductors=cable_geometry_list(cable_geometry)%n_conductors
      
      conductor_count=conductor_count+n_conductors
      filter_count=filter_count+cable_geometry_list(cable_geometry)%n_filters
      
    end do ! next cable
    
! Allocate memory for L, C and R matrices

    total_n_bundle_conductors=conductor_count
    bundle_segment_geometry_list(bundle_segment_geometry)%n_conductors=total_n_bundle_conductors
    bundle_segment_geometry_list(bundle_segment_geometry)%n_filters=filter_count
        
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%L(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%C(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%R(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Zlink(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Ylink(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%ZLStub(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Yf(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Tv(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Ti(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%SC(1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%excitation_function(1:conductor_count) )
    
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%filter_number(1:conductor_count,1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Sfilter(1:filter_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Zfilter(1:filter_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%Z_f(1:filter_count) )

! set all matrices to zero
    bundle_segment_geometry_list(bundle_segment_geometry)%L(1:conductor_count,1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%C(1:conductor_count,1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%R(1:conductor_count,1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%Zlink(1:conductor_count,1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%Ylink(1:conductor_count,1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%ZLstub(1:conductor_count,1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%Yf(1:conductor_count,1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%Tv(1:conductor_count,1:conductor_count)=0
    bundle_segment_geometry_list(bundle_segment_geometry)%Ti(1:conductor_count,1:conductor_count)=0
    bundle_segment_geometry_list(bundle_segment_geometry)%SC(1:conductor_count)=0
    bundle_segment_geometry_list(bundle_segment_geometry)%excitation_function(1:conductor_count)=0
    
    bundle_segment_geometry_list(bundle_segment_geometry)%filter_number(1:conductor_count,1:conductor_count)=0
   
! Allocate memory for segment geometry
        
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%xc(1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%yc(1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%rc(1:conductor_count) )
    ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry)%ri(1:conductor_count) )
        
    bundle_segment_geometry_list(bundle_segment_geometry)%xc(1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%yc(1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%rc(1:conductor_count)=0d0
    bundle_segment_geometry_list(bundle_segment_geometry)%ri(1:conductor_count)=0d0

! set the shielded cable part of the L and C matrices    
    
    conductor_count=0
    filter_count=0
    do cable_loop=1,bundle_segment_geometry_list(bundle_segment_geometry)%n_cables

      cable=bundle_segment_geometry_list(bundle_segment_geometry)%cable_list(cable_loop)
      cable_geometry=cable_list(cable)%cable_geometry_number
      n_conductors=cable_geometry_list(cable_geometry)%n_conductors
      n_shielded_conductors=cable_geometry_list(cable_geometry)%n_shielded_conductors
      n_external_conductors=cable_geometry_list(cable_geometry)%n_external_conductors
      
      n_filters=cable_geometry_list(cable_geometry)%n_filters
      
! Copy local cable internal L, C, Tv, Ti and SC data 
! into bundle_segment_geometry_list(bundle_segment_geometry)%L, C, Tv, Ti and Sc data

      do cable_row=1,n_conductors
      
        bundle_row=conductor_count+cable_row   
	   
        do cable_col=1,n_conductors
	
          bundle_col=conductor_count+cable_col   
          
	  bundle_segment_geometry_list(bundle_segment_geometry)%L(bundle_row,bundle_col)=	&
	         cable_geometry_list(cable_geometry)%L_internal(cable_row,cable_col)
		 
	  bundle_segment_geometry_list(bundle_segment_geometry)%C(bundle_row,bundle_col)=	&
	         cable_geometry_list(cable_geometry)%C_internal(cable_row,cable_col)
		 
	  bundle_segment_geometry_list(bundle_segment_geometry)%R(bundle_row,bundle_col)=	&
	         cable_geometry_list(cable_geometry)%R_internal(cable_row,cable_col)
		 
	  bundle_segment_geometry_list(bundle_segment_geometry)%Tv(bundle_row,bundle_col)=	&
	         cable_geometry_list(cable_geometry)%Tv(cable_row,cable_col)
		 
	  bundle_segment_geometry_list(bundle_segment_geometry)%Ti(bundle_row,bundle_col)=	&
	         cable_geometry_list(cable_geometry)%Ti(cable_row,cable_col)
		 
	  bundle_segment_geometry_list(bundle_segment_geometry)%filter_number(bundle_row,bundle_col)=	&
	         cable_geometry_list(cable_geometry)%filter_number(cable_row,cable_col)+filter_count
	  
        end do ! next col
		 
	bundle_segment_geometry_list(bundle_segment_geometry)%SC(bundle_row)=	&
	       cable_geometry_list(cable_geometry)%SC(cable_row)
	
      end do ! next row
      
      conductor_count=conductor_count+n_conductors

! copy filters from the cable_geometry_list to the bundle_segment_geometry_list
      do filter=1,n_filters
        bundle_segment_geometry_list(bundle_segment_geometry)%Sfilter(filter+filter_count)=	&
	                  cable_geometry_list(cable_geometry)%Sfilter(filter)
        bundle_segment_geometry_list(bundle_segment_geometry)%Zfilter(filter+filter_count)=	&
	                  cable_geometry_list(cable_geometry)%Zfilter(filter)
        bundle_segment_geometry_list(bundle_segment_geometry)%Z_f(filter+filter_count)=    &
	                  cable_geometry_list(cable_geometry)%Z_f(filter)
      end do
      filter_count=filter_count+n_filters

    end do ! next cable

! set up the structure for LC matrix calculation, one element for each cable 
! for the initial process of determining the cable arrangement in the bundle

    n_LC_conductors=bundle_segment_geometry_list(bundle_segment_geometry)%n_cables
        
    ALLOCATE( cable_radius(1:n_LC_conductors+1) ) ! include extra element for TLM return conductor radius
    ALLOCATE( cable_xc(1:n_LC_conductors+1) )     ! include extra element for TLM return conductor
    ALLOCATE( cable_yc(1:n_LC_conductors+1) )     ! include extra element for TLM return conductor
    
    cable_radius(1:n_LC_conductors+1)=0d0
    cable_xc(1:n_LC_conductors+1)=0d0
    cable_yc(1:n_LC_conductors+1)=0d0
    
    do cable_loop=1,bundle_segment_geometry_list(bundle_segment_geometry)%n_cables
    
      cable=bundle_segment_geometry_list(bundle_segment_geometry)%cable_list(cable_loop)
      cable_geometry=cable_list(cable)%cable_geometry_number
      cable_radius(cable_loop)=cable_geometry_list(cable_geometry)%cable_offset_radius
      
    end do
    
    write(cable_info_file_unit,*)'CALL create_cable_geometry'
    write(cable_info_file_unit,*)'matrix dimension',n_LC_conductors+1
    flush(cable_info_file_unit) 
    
    CALL create_bundle_cross_section_geometry(	n_LC_conductors,cable_radius,	&
      						cable_xc,cable_yc)

    
    write(cable_info_file_unit,*)''
    write(cable_info_file_unit,*)' cable       cable        cable         cable    '
    write(cable_info_file_unit,*)' number     radius	     xc	           yc	   '
    do row=1,n_LC_conductors+1
      write(cable_info_file_unit,8000)row,cable_radius(row),cable_xc(row),cable_yc(row)
    end do
     
! We now have the arrangement of the cables in the bundle, now calculate the
! arrangement of the individual conductors in the bundle
    
! set up the structure for LC matrix calculation, each cable contributes n_external conductors now

    n_LC_conductors=0
    
    do cable_loop=1,bundle_segment_geometry_list(bundle_segment_geometry)%n_cables

      cable=bundle_segment_geometry_list(bundle_segment_geometry)%cable_list(cable_loop)
      cable_geometry=cable_list(cable)%cable_geometry_number
      
      n_external_conductors=cable_geometry_list(cable_geometry)%n_external_conductors
      
      n_LC_conductors=n_LC_conductors+n_external_conductors
      
    end do
    
    write(cable_info_file_unit,*)'number of external conductors=',n_LC_conductors
    
    write(cable_info_file_unit,*)'CALL create_conductor_section_geometry'
    write(cable_info_file_unit,*)'matrix dimension',n_LC_conductors+1
    flush(cable_info_file_unit) 
    
    ALLOCATE( conductor_radius(1:n_LC_conductors+1) ) ! include extra element for TLM return conductor
    ALLOCATE( conductor_xc(1:n_LC_conductors+1) )     ! include extra element for TLM return conductor
    ALLOCATE( conductor_yc(1:n_LC_conductors+1) )     ! include extra element for TLM return conductor
    ALLOCATE( dielectric_radius(1:n_LC_conductors+1) )
    ALLOCATE( dielectric_permittivity(1:n_LC_conductors+1) )
    ALLOCATE( local_L(1:n_LC_conductors,1:n_LC_conductors) )
    ALLOCATE( local_C(1:n_LC_conductors,1:n_LC_conductors) )
    ALLOCATE( local_LC_conductor_to_segment_conductor(1:n_LC_conductors) )
    
    conductor_xc(1:n_LC_conductors+1)=0d0
    conductor_yc(1:n_LC_conductors+1)=0d0
    dielectric_radius(1:n_LC_conductors+1)=0d0
    dielectric_permittivity(1:n_LC_conductors+1)=0d0
    local_L(1:n_LC_conductors,1:n_LC_conductors)=0d0
    local_C(1:n_LC_conductors,1:n_LC_conductors)=0d0
    local_LC_conductor_to_segment_conductor(1:n_LC_conductors)=0
    
    conductor_count=0
    total_conductor_count=0
    do cable_loop=1,bundle_segment_geometry_list(bundle_segment_geometry)%n_cables
    
      cable=bundle_segment_geometry_list(bundle_segment_geometry)%cable_list(cable_loop)
      cable_geometry=cable_list(cable)%cable_geometry_number
      
      n_external_conductors=cable_geometry_list(cable_geometry)%n_external_conductors
      n_shielded_conductors=cable_geometry_list(cable_geometry)%n_shielded_conductors
	
      cable_centre_x=cable_xc(cable_loop)
      cable_centre_y=cable_yc(cable_loop)
      
      do i=1,n_external_conductors

! external conductors for L and C matrix calculation      
        conductor_count=conductor_count+1	
        conductor_radius(conductor_count)=cable_geometry_list(cable_geometry)%external_conductor_radius(i)
        dielectric_radius(conductor_count)=cable_geometry_list(cable_geometry)%external_dielectric_radius(i)	
        conductor_xc(conductor_count)=cable_geometry_list(cable_geometry)%external_conductor_xc(i)+cable_xc(cable_loop)
        conductor_yc(conductor_count)=cable_geometry_list(cable_geometry)%external_conductor_yc(i)+cable_yc(cable_loop)	
        dielectric_permittivity(conductor_count)=cable_geometry_list(cable_geometry)%external_dielectric_permittivity(i)

! external conductors for bundle_segment geometry specification    
        total_conductor_count=total_conductor_count+1
		
        bundle_segment_geometry_list(bundle_segment_geometry)%xc(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%external_conductor_xc(i)+cable_xc(cable_loop)
        bundle_segment_geometry_list(bundle_segment_geometry)%yc(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%external_conductor_yc(i)+cable_yc(cable_loop)
        bundle_segment_geometry_list(bundle_segment_geometry)%rc(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%external_conductor_radius(i)
        bundle_segment_geometry_list(bundle_segment_geometry)%ri(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%external_dielectric_radius(i)
	
      end do ! next external conductor in this cable
      
      do i=1,n_shielded_conductors

! internal conductors for bundle_segment geometry specification    
        total_conductor_count=total_conductor_count+1
		
        bundle_segment_geometry_list(bundle_segment_geometry)%xc(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%shielded_conductor_xc(i)+cable_xc(cable_loop)
        bundle_segment_geometry_list(bundle_segment_geometry)%yc(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%shielded_conductor_yc(i)+cable_yc(cable_loop)
        bundle_segment_geometry_list(bundle_segment_geometry)%rc(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%shielded_conductor_radius(i)
        bundle_segment_geometry_list(bundle_segment_geometry)%ri(total_conductor_count)=	&
	                            cable_geometry_list(cable_geometry)%shielded_dielectric_radius(i)
	
      end do ! next internal conductor in this cable
      
    end do

! radius of cable bundle in the cell    
    conductor_count=n_LC_conductors+1
    conductor_radius(conductor_count)=cable_radius(bundle_segment_geometry_list(bundle_segment_geometry)%n_cables+1)
    dielectric_radius(conductor_count)=conductor_radius(conductor_count)

    write(cable_info_file_unit,*)'Radius of cable bundle=',conductor_radius(conductor_count)

    bundle_segment_geometry_list(bundle_segment_geometry)%cable_bundle_radius=conductor_radius(conductor_count)
    
    write(cable_info_file_unit,*)''
    write(cable_info_file_unit,*)'conductor  conductor    conductor     conductor    dielectric     dielectric'
    write(cable_info_file_unit,*)' number     radius	     xc	           yc	      radius      permittivity'
    do row=1,n_LC_conductors+1
      write(cable_info_file_unit,8000)row,conductor_radius(row),conductor_xc(row),conductor_yc(row),	&
                              dielectric_radius(row),dielectric_permittivity(row)
    end do
8000 format(I7,5E14.4)
    
    flush(cable_info_file_unit)
    
    DEALLOCATE( cable_radius )
    DEALLOCATE( cable_xc )
    DEALLOCATE( cable_yc )
     
! Call PUL_LC_CALC through the interface subroutine PUL_LC_TLM_INTERFACE
! This returns the per-unit length Inductance and Capacitance matrices
  
    CALL PUL_LC_tlm_interface(n_LC_conductors,conductor_radius,conductor_xc,conductor_yc,	&
                              dielectric_radius,dielectric_permittivity,local_L,local_C,	&
			      reference_radius,reference_radius_L,reference_radius_C)
    
    write(cable_info_file_unit,*)'L_external'    
    do row=1,n_LC_conductors
      write(cable_info_file_unit,8010)(local_L(row,col),col=1,n_LC_conductors)
    end do ! next row
    
    write(cable_info_file_unit,*)'C_external'    
    do row=1,n_LC_conductors
      write(cable_info_file_unit,8010)(local_C(row,col),col=1,n_LC_conductors)
    end do ! next row
    
    write(cable_info_file_unit,*)''
    
8010  format(100E14.4)

    flush(cable_info_file_unit)
			      
! work out the mapping from the local LC matrices to the bundle segment LC matrices

    conductor_count=0
    external_conductor_count=0
    do cable_loop=1,bundle_segment_geometry_list(bundle_segment_geometry)%n_cables
    
      cable=bundle_segment_geometry_list(bundle_segment_geometry)%cable_list(cable_loop)
      cable_geometry=cable_list(cable)%cable_geometry_number
    
! the unshielded conductors take the last conductor numbers for the cable            
      n_conductors=cable_geometry_list(cable_geometry)%n_conductors
      n_shielded_conductors=cable_geometry_list(cable_geometry)%n_shielded_conductors  
      n_external_conductors=cable_geometry_list(cable_geometry)%n_external_conductors
      
      do i=1,n_shielded_conductors
        conductor_count=conductor_count+1
      end do
      
      do i=1,n_external_conductors
        conductor_count=conductor_count+1
	external_conductor_count=external_conductor_count+1
        local_LC_conductor_to_segment_conductor(external_conductor_count)=conductor_count
      end do
      
    end do
    
! Put the unshielded cable L and C matrix elements into the bundle_segment_geometry L and C matrices
! and scale to give per segment L and C ie for a segment of length dl/2  
! The in-cell inductance is removed from the cable inductance diagonal terms here if this was requested
! in the input file
 
    segment_length=dl/2d0

    do row=1,n_LC_conductors
      bundle_row=local_LC_conductor_to_segment_conductor(row)
      do col=1,n_LC_conductors
      
        bundle_col=local_LC_conductor_to_segment_conductor(col)
        bundle_segment_geometry_list(bundle_segment_geometry)%L(bundle_row,bundle_col)=	&
	                                                local_L(row,col)*segment_length
        bundle_segment_geometry_list(bundle_segment_geometry)%C(bundle_row,bundle_col)=	&
	                                                local_C(row,col)*segment_length

      end do
    end do

! The in-cell inductance is removed from the cable inductance diagonal terms here if this was requested
! in the input file
! note that the inductance is only corrected for unshielded conductors
    if (Cable_LC_Correction_type.EQ.LC_correction_type_subtract_cell_inductance) then

      do row=1,n_LC_conductors
        bundle_row=local_LC_conductor_to_segment_conductor(row)      
        bundle_col=bundle_row
	
        if (bundle_segment_geometry_list(bundle_segment_geometry)%SC(bundle_row).EQ.0) then
	
          bundle_segment_geometry_list(bundle_segment_geometry)%L(bundle_row,bundle_col)=	&
                bundle_segment_geometry_list(bundle_segment_geometry)%L(bundle_row,bundle_col)	&
	        -(Z0*dt/16d0)
		
        end if
	
      end do
    
    end if  ! subtract in-cell inductance
    
    bundle_segment_geometry_list(bundle_segment_geometry)%TLM_cell_equivalent_radius=1.08D0*dl/2d0
    bundle_segment_geometry_list(bundle_segment_geometry)%TLM_reference_radius_rL=reference_radius_L
    bundle_segment_geometry_list(bundle_segment_geometry)%TLM_reference_radius_rC=reference_radius_C

    write(cable_info_file_unit,*)'Cable bundle radius         =',	&
          bundle_segment_geometry_list(bundle_segment_geometry)%cable_bundle_radius
    write(cable_info_file_unit,*)'Half TLM cell size          =',dl/2d0
    write(cable_info_file_unit,*)'Equivalent cylider radius   =',1.08D0*dl/2d0
    write(cable_info_file_unit,*)'Inductance reference radius =',	&
          bundle_segment_geometry_list(bundle_segment_geometry)%TLM_reference_radius_rL
    write(cable_info_file_unit,*)'Capacitance reference radius=',	&
          bundle_segment_geometry_list(bundle_segment_geometry)%TLM_reference_radius_rC

! Check that the cable matrices describe a physical system    
    
    CALL LC_checks(bundle_segment_geometry_list(bundle_segment_geometry)%L,	&
                   bundle_segment_geometry_list(bundle_segment_geometry)%C,	&
		   segment_length,total_n_bundle_conductors,total_n_bundle_conductors)

! Calculate the TLM cable update matrices

    CALL calculate_TLM_cable_matrices(bundle_segment_geometry_list(bundle_segment_geometry)%L,	&
                                      bundle_segment_geometry_list(bundle_segment_geometry)%C,	&
                                      bundle_segment_geometry_list(bundle_segment_geometry)%Zlink,	&
                                      bundle_segment_geometry_list(bundle_segment_geometry)%Zlstub,	&
				      total_n_bundle_conductors) ! was n_LC_conductors... IN ERROR

! Check the TLM cable update matrices for stability and correct if required		   
    
    write_to_cable_info_file_flag=.TRUE.
    CALL TLM_Cable_checks(bundle_segment_geometry_list(bundle_segment_geometry)%Zlink,	&
                          bundle_segment_geometry_list(bundle_segment_geometry)%ZLstub,	&
			  total_n_bundle_conductors,write_to_cable_info_file_flag) ! was n_LC_conductors... IN ERROR
			  

    call dsvd_invert(bundle_segment_geometry_list(bundle_segment_geometry)%Zlink	&
                     ,total_n_bundle_conductors,total_n_bundle_conductors,	&
		     bundle_segment_geometry_list(bundle_segment_geometry)%Ylink,total_n_bundle_conductors)

! This Zf, Yf calculation is repeated in cable_bundles.F90  on a segment basis 
! It is done here so that the matrices can be looked at in GGI_TLM_cable_model_checks
! We may rationalise this as it is wasteful.   
    ALLOCATE( Zf(1:total_n_bundle_conductors,1:total_n_bundle_conductors) )
        
    Zf(:,:)=bundle_segment_geometry_list(bundle_segment_geometry)%Zlink(:,:)+	&
            bundle_segment_geometry_list(bundle_segment_geometry)%ZLstub(:,:)+	&
	    bundle_segment_geometry_list(bundle_segment_geometry)%R(:,:)

! include fast impedance from filter function into Zf

    do row=1,total_n_bundle_conductors
      do col=1,total_n_bundle_conductors
      
        filter=bundle_segment_geometry_list(bundle_segment_geometry)%filter_number(row,col)
	if (filter.ne.0) then ! there is a filter on this matrix element
	
	  Zf(row,col)=Zf(row,col)+	&
	    bundle_segment_geometry_list(bundle_segment_geometry)%Z_f(filter) 
			
	end if
      end do
    end do
    		
    call dsvd_invert(Zf,total_n_bundle_conductors,total_n_bundle_conductors,	&
                    bundle_segment_geometry_list(bundle_segment_geometry)%Yf,total_n_bundle_conductors) 
			  
    DEALLOCATE( conductor_radius )
    DEALLOCATE( conductor_xc )
    DEALLOCATE( conductor_yc )
    DEALLOCATE( dielectric_radius )
    DEALLOCATE( dielectric_permittivity )
    DEALLOCATE( local_L )
    DEALLOCATE( local_C )     
    DEALLOCATE( Zf )     
    DEALLOCATE( local_LC_conductor_to_segment_conductor )
     
  end do ! next bundle segment geometry

! Distribute the L and C matrices to the required bundle segments

  write(cable_info_file_unit,*)'       Bundle      Bundle     n_conductors'
  write(cable_info_file_unit,*)'      segment     geometry'
  
  do bundle_segment=1,n_bundle_segments
    
    bundle_segment_geometry=bundle_segment_list(bundle_segment)%bundle_segment_geometry
	
! Allocate memory for L, C and R matrices
    bundle_segment_list(bundle_segment)%n_conductors=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%n_conductors
    bundle_segment_list(bundle_segment)%n_filters=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%n_filters
    
    n_conductors=bundle_segment_list(bundle_segment)%n_conductors
    n_filters=bundle_segment_list(bundle_segment)%n_filters
    
    write(cable_info_file_unit,*)bundle_segment,bundle_segment_geometry,n_conductors
    
    ALLOCATE( bundle_segment_list(bundle_segment)%L(1:n_conductors,1:n_conductors) )
    ALLOCATE( bundle_segment_list(bundle_segment)%C(1:n_conductors,1:n_conductors) )
    ALLOCATE( bundle_segment_list(bundle_segment)%R(1:n_conductors,1:n_conductors) )
    ALLOCATE( bundle_segment_list(bundle_segment)%Tv(1:n_conductors,1:n_conductors) )
    ALLOCATE( bundle_segment_list(bundle_segment)%Ti(1:n_conductors,1:n_conductors) )
    ALLOCATE( bundle_segment_list(bundle_segment)%SC(1:n_conductors) )
    ALLOCATE( bundle_segment_list(bundle_segment)%excitation_function(1:n_conductors) )
    
    ALLOCATE( bundle_segment_list(bundle_segment)%filter_number(1:n_conductors,1:n_conductors) )
    ALLOCATE( bundle_segment_list(bundle_segment)%Sfilter(1:n_filters) )
    ALLOCATE( bundle_segment_list(bundle_segment)%Zfilter(1:n_filters) )
    ALLOCATE( bundle_segment_list(bundle_segment)%Z_f(1:n_filters) )
    
    ALLOCATE( bundle_segment_list(bundle_segment)%direction_sign_list(1:n_conductors) )
    bundle_segment_list(bundle_segment)%direction_sign_list(1:n_conductors)=0

! copy matrices from bundle_segment_geometry_list to bundle_segment_list
    bundle_segment_list(bundle_segment)%L(1:n_conductors,1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%L(1:n_conductors,1:n_conductors)
    
    bundle_segment_list(bundle_segment)%C(1:n_conductors,1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%C(1:n_conductors,1:n_conductors)
    
    bundle_segment_list(bundle_segment)%R(1:n_conductors,1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%R(1:n_conductors,1:n_conductors)
    
    bundle_segment_list(bundle_segment)%Tv(1:n_conductors,1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%Tv(1:n_conductors,1:n_conductors)
    
    bundle_segment_list(bundle_segment)%Ti(1:n_conductors,1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%Ti(1:n_conductors,1:n_conductors)
    
    bundle_segment_list(bundle_segment)%SC(1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%SC(1:n_conductors)
    
    bundle_segment_list(bundle_segment)%excitation_function(1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%excitation_function(1:n_conductors)

    
    bundle_segment_list(bundle_segment)%filter_number(1:n_conductors,1:n_conductors)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%filter_number(1:n_conductors,1:n_conductors)
    
    bundle_segment_list(bundle_segment)%Sfilter(1:n_filters)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%Sfilter(1:n_filters)
    
    bundle_segment_list(bundle_segment)%Zfilter(1:n_filters)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%Zfilter(1:n_filters)
    
    bundle_segment_list(bundle_segment)%Z_f(1:n_filters)=	&
    bundle_segment_geometry_list(bundle_segment_geometry)%Z_f(1:n_filters)

  end do ! next bundle segment

! set direction_sign_list on a conductor basis in the bundle segment

  ALLOCATE( segment_conductor_count(1:n_bundle_segments) )
  segment_conductor_count(1:n_bundle_segments)=0
  
  do cable=1,n_cables
  
    cable_geometry=cable_list(cable)%cable_geometry_number
    n_conductors=cable_geometry_list(cable_geometry)%n_conductors
  
    do cable_segment=1,cable_list(cable)%number_of_cable_segments
    
      direction_sign=cable_list(cable)%direction_sign_list(cable_segment)
      bundle_segment=cable_list(cable)%bundle_segment_list(cable_segment)
      
      first_conductor=segment_conductor_count(bundle_segment)+1
      last_conductor=segment_conductor_count(bundle_segment)+n_conductors
      
      bundle_segment_list(bundle_segment)%direction_sign_list(first_conductor:last_conductor)=direction_sign

      segment_conductor_count(bundle_segment)=last_conductor
      
    end do ! next segment
    
  end do ! next cable
  
! check - the conductor count should be equal to n_conductors for all the bundle_segments now if all
!         direction signs have been correctly allocated on a conductor basis.
	  
  do bundle_segment=1,n_bundle_segments
    
    if (segment_conductor_count(bundle_segment).ne.bundle_segment_list(bundle_segment)%n_conductors) then
      write(*,*)'Error setting direction_sign_list on a conductor basis in the bundle_segment_list'
      write(*,*)'Conductor counting discrepancy on bundle segment',bundle_segment
      write(*,*)'segment_conductor_count(bundle_segment)=',segment_conductor_count(bundle_segment)
      write(*,*)'bundle_segment_list(bundle_segment)%n_conductors=',bundle_segment_list(bundle_segment)%n_conductors
      STOP
    end if
    
  end do  

  DEALLOCATE( segment_conductor_count )
  
  CALL write_line('FINISHED: create_bundle_LCRG_matrices',0,output_to_screen_flag)
    
  RETURN
  
  
END SUBROUTINE create_bundle_LCRG_matrices
