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
!SUBROUTINE read_cable_model
!
! NAME
!     SUBROUTINE read_cable_model
!
! DESCRIPTION
!      read the cable model data from file
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!
!
SUBROUTINE read_cable_model(read_data_for_computation_only)

USE TLM_general
USE Cables
USE File_information
USE filter_types
USE filter_functions

IMPLICIT NONE

! variables passed to subroutine

logical :: read_data_for_computation_only

! local variables

  integer 	:: segment
  integer 	:: segment_geometry
  integer	:: cell,cell_face
  integer	:: cable
  integer	:: cable_geometry
  integer	:: output
  integer	:: conductor
  integer	:: n_conductors
  integer	:: n_filters
  integer	:: n_external,n_internal,n_shielded
  integer	:: face
  integer	:: row,col
  integer	:: n_rows,n_cols
  integer	:: i
  integer	:: filter
  
  character*256	:: cable_model_file_name

  logical	:: file_exists

! START

  CALL write_line('CALLED: read_cable_model',0,output_to_screen_flag)

! check for the existance of a cable_model file

  cable_model_file_name=trim(problem_name)//cable_model_file_extn

  inquire(file=trim(cable_model_file_name),exist=file_exists)
  
  if (.NOT.file_exists) then

! set cable bundle numbers to zero and return  
    CALL write_line('Cable_model file not found',0,output_to_screen_flag)
    n_bundle_segments=0
    n_cell_centre_junctions=0
    n_face_junctions=0
    RETURN
    
  end if
  
! Open cable model file

  CALL open_file(cable_model_file_unit,cable_model_file_extn)
  
! READ BUNDLE_SEGMENT_LIST
  CALL write_line('Read bundle_segment_list',0,output_to_screen_flag)
  
  read(cable_model_file_unit,*) ! Read comment line
  
  read(cable_model_file_unit,*)n_bundle_segments
  
  ALLOCATE( bundle_segment_list(1:n_bundle_segments) )
  
  do segment=1,n_bundle_segments
    
    read(cable_model_file_unit,*)bundle_segment_list(segment)%cable_segment%segment_point(1)%cell%i,	&
    				 bundle_segment_list(segment)%cable_segment%segment_point(1)%cell%j,   &
    				 bundle_segment_list(segment)%cable_segment%segment_point(1)%cell%k,   &
    				 bundle_segment_list(segment)%cable_segment%segment_point(1)%point,    &
    				 bundle_segment_list(segment)%cable_segment%segment_point(2)%cell%i,   &
    				 bundle_segment_list(segment)%cable_segment%segment_point(2)%cell%j,   &
    				 bundle_segment_list(segment)%cable_segment%segment_point(2)%cell%k,   &
    				 bundle_segment_list(segment)%cable_segment%segment_point(2)%point
    
    read(cable_model_file_unit,*)bundle_segment_list(segment)%n_cables
    
    ALLOCATE( bundle_segment_list(segment)%cable_list(1:bundle_segment_list(segment)%n_cables) )
    
    do cable=1,bundle_segment_list(segment)%n_cables    
      read(cable_model_file_unit,*)bundle_segment_list(segment)%cable_list(cable)   
    end do ! next cable
    
    read(cable_model_file_unit,*)bundle_segment_list(segment)%n_conductors   
    read(cable_model_file_unit,*)bundle_segment_list(segment)%n_filters
    
    n_conductors=bundle_segment_list(segment)%n_conductors 
    n_filters=bundle_segment_list(segment)%n_filters
    
    ALLOCATE( bundle_segment_list(segment)%direction_sign_list(1:n_conductors) )
    
    do row=1,n_conductors
      read(cable_model_file_unit,*)bundle_segment_list(segment)%direction_sign_list(row)
    end do ! next conductor
    
    read(cable_model_file_unit,*)bundle_segment_list(segment)%bundle_segment_geometry
    
    ALLOCATE( bundle_segment_list(segment)%L(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%C(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%R(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%Tv(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%Ti(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%SC(1:n_conductors) )
    ALLOCATE( bundle_segment_list(segment)%excitation_function(1:n_conductors) )
    
    ALLOCATE( bundle_segment_list(segment)%filter_number(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_list(segment)%Sfilter(1:n_filters) )
    ALLOCATE( bundle_segment_list(segment)%Zfilter(1:n_filters) )
    ALLOCATE( bundle_segment_list(segment)%Z_f(1:n_filters) )
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_list(segment)%L(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_list(segment)%C(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_list(segment)%R(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_list(segment)%Tv(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_list(segment)%Ti(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      read(cable_model_file_unit,*)bundle_segment_list(segment)%SC(row)           
    end do ! next row
    
    do row=1,n_conductors
      read(cable_model_file_unit,*)bundle_segment_list(segment)%excitation_function(row)           
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_list(segment)%filter_number(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_filters
      read(cable_model_file_unit,*)
      CALL read_Sfilter(bundle_segment_list(segment)%Sfilter(row),cable_model_file_unit)        
    end do ! next row
    
    do row=1,n_filters
      read(cable_model_file_unit,*)          
      read(cable_model_file_unit,*)bundle_segment_list(segment)%Z_f(row)
      read(cable_model_file_unit,*)        
      CALL read_Zfilter(bundle_segment_list(segment)%Zfilter(row),cable_model_file_unit)
    end do ! next row
	
  end do ! next bundle segment
  
! READ BUNDLE_SEGMENT_GEOMETYRY LIST
  
  read(cable_model_file_unit,*) ! Read comment line

  CALL write_line('Read bundle_segment_geometry_list',0,output_to_screen_flag)

  read(cable_model_file_unit,*)n_bundle_segment_geometries
  ALLOCATE( bundle_segment_geometry_list(1:n_bundle_segment_geometries) )

  do segment_geometry=1,n_bundle_segment_geometries
        
    read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%n_cables
    
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%cable_list(1:bundle_segment_geometry_list(segment_geometry)%n_cables) )
    
    do cable=1,bundle_segment_geometry_list(segment_geometry)%n_cables    
      read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%cable_list(cable)   
    end do ! next cable
    
    read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%n_conductors

    n_conductors=bundle_segment_geometry_list(segment_geometry)%n_conductors 

    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%xc(1:n_conductors) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%yc(1:n_conductors) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%rc(1:n_conductors) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%ri(1:n_conductors) )
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%xc(row),	&
                                   bundle_segment_geometry_list(segment_geometry)%yc(row),	&
                                   bundle_segment_geometry_list(segment_geometry)%rc(row),	&
                                   bundle_segment_geometry_list(segment_geometry)%ri(row)    
    end do ! next row   
    
    read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%cable_bundle_radius
    read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%TLM_cell_equivalent_radius
    read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%TLM_reference_radius_rL
    read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%TLM_reference_radius_rC
    
    read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%n_filters
    
    n_filters=bundle_segment_geometry_list(segment_geometry)%n_filters
    
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%L(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%C(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%R(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Zlink(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Ylink(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%ZLstub(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Yf(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Tv(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Ti(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%SC(1:n_conductors) )
    
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%filter_number(1:n_conductors,1:n_conductors ) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Sfilter(1:n_filters) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Zfilter(1:n_filters) )
    ALLOCATE( bundle_segment_geometry_list(segment_geometry)%Z_f(1:n_filters) )
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%L(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%C(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%R(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Zlink(row,col)         
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Ylink(row,col)         
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%ZLstub(row,col)         
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Yf(row,col)         
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Tv(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Ti(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_conductors
      read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%SC(row)           
    end do ! next row
    
    do row=1,n_conductors
      do col=1,n_conductors
        read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%filter_number(row,col)           
      end do ! next col
    end do ! next row
    
    do row=1,n_filters
      read(cable_model_file_unit,*)
      CALL read_Sfilter(bundle_segment_geometry_list(segment_geometry)%Sfilter(row),cable_model_file_unit)        
    end do ! next row
    
    do row=1,n_filters
      read(cable_model_file_unit,*)          
      read(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Z_f(row)
      read(cable_model_file_unit,*)        
      CALL read_Zfilter(bundle_segment_geometry_list(segment_geometry)%Zfilter(row),cable_model_file_unit)
    end do ! next row
 	
  end do ! next bundle segment_geometry

! READ CELL_CENTRE_JUNCTION_LIST
  CALL write_line('Read cell centre junction list',0,output_to_screen_flag)
  
  read(cable_model_file_unit,*) ! Read comment line

  read(cable_model_file_unit,*)n_cell_centre_junctions

  ALLOCATE( cell_centre_junction_list(1:n_cell_centre_junctions) )

  do cell=1,n_cell_centre_junctions
  
    read(cable_model_file_unit,*)cell_centre_junction_list(cell)%cell_point%cell%i,	&
                                 cell_centre_junction_list(cell)%cell_point%cell%j, 	&
                                 cell_centre_junction_list(cell)%cell_point%cell%k, 	&
                                 cell_centre_junction_list(cell)%cell_point%point
				     
    read(cable_model_file_unit,*)cell_centre_junction_list(cell)%n_internal_connection_nodes
    read(cable_model_file_unit,*)cell_centre_junction_list(cell)%n_segments
    
    if (cell_centre_junction_list(cell)%n_segments.NE.6) GOTO 9000   
    ALLOCATE( cell_centre_junction_list(cell)%segment_list(1:6) )
    ALLOCATE( cell_centre_junction_list(cell)%n_external_conductors(1:7) )
    ALLOCATE( cell_centre_junction_list(cell)%P_matrix_list(1:7) )
    
    read(cable_model_file_unit,*)cell_centre_junction_list(cell)%segment_list(1:6)
    read(cable_model_file_unit,*)cell_centre_junction_list(cell)%n_external_conductors(1:7)

    read(cable_model_file_unit,*)cell_centre_junction_list(cell)%number_of_cable_junctions
    
    ALLOCATE( cell_centre_junction_list(cell)%cable_junction_list(1:cell_centre_junction_list(cell)%number_of_cable_junctions) )

    do i=1,cell_centre_junction_list(cell)%number_of_cable_junctions
      read(cable_model_file_unit,*)cell_centre_junction_list(cell)%cable_junction_list(i)
    end do
    
    n_internal=cell_centre_junction_list(cell)%n_internal_connection_nodes
    do face=1,7
    
      n_external=cell_centre_junction_list(cell)%n_external_conductors(face)
      if (n_external.ne.0) then
	
        ALLOCATE( cell_centre_junction_list(cell)%P_matrix_list(face)%P(1:n_internal,1:n_external) )
        do row=1,cell_centre_junction_list(cell)%n_internal_connection_nodes
	
          read(cable_model_file_unit,*)cell_centre_junction_list(cell)%P_matrix_list(face)%P(row,1:n_external)
	  
        end do
	
      end if ! n_external.gt.0
      
    end do ! next face    
    
! Internal impedance stuff    
    read(cable_model_file_unit,*)cell_centre_junction_list(cell)%n_internal_impedance_filters
    n_filters=cell_centre_junction_list(cell)%n_internal_impedance_filters
    
    ALLOCATE( cell_centre_junction_list(cell)%Sfilter(1:n_filters) )
    ALLOCATE( cell_centre_junction_list(cell)%Zfilter(1:n_filters) )
    ALLOCATE( cell_centre_junction_list(cell)%Z_f(1:n_filters) )
   
    do filter=1,cell_centre_junction_list(cell)%n_internal_impedance_filters

      read(cable_model_file_unit,*)  
      CALL read_Sfilter(cell_centre_junction_list(cell)%Sfilter(filter),cable_model_file_unit)
      read(cable_model_file_unit,*)        
      read(cable_model_file_unit,*)cell_centre_junction_list(cell)%Z_f(filter)
      read(cable_model_file_unit,*)   
      CALL read_Zfilter(cell_centre_junction_list(cell)%Zfilter(filter),cable_model_file_unit)
      
    end do ! next filter
      
  end do ! next cable cell

! READ FACE_JUNCTION_LIST
  CALL write_line('Read face junction list',0,output_to_screen_flag)
  
  read(cable_model_file_unit,*) ! Read comment line
  
  read(cable_model_file_unit,*)n_face_junctions

  ALLOCATE( face_junction_list(1:n_face_junctions) )

  do cell_face=1,n_face_junctions
  
    read(cable_model_file_unit,*)face_junction_list(cell_face)%cell_point%cell%i,	&
                                 face_junction_list(cell_face)%cell_point%cell%j, 	&
                                 face_junction_list(cell_face)%cell_point%cell%k, 	&
                                 face_junction_list(cell_face)%cell_point%point
				     
    read(cable_model_file_unit,*)face_junction_list(cell_face)%n_internal_connection_nodes
    
    ALLOCATE( face_junction_list(cell_face)%BC(1:face_junction_list(cell_face)%n_internal_connection_nodes) )
    
    do row=1,face_junction_list(cell_face)%n_internal_connection_nodes
      read(cable_model_file_unit,*)face_junction_list(cell_face)%BC(row)
    end do
    
    read(cable_model_file_unit,*)face_junction_list(cell_face)%n_segments
    
    if (face_junction_list(cell_face)%n_segments.NE.2) GOTO 9010
    ALLOCATE( face_junction_list(cell_face)%segment_list(1:2) )
    ALLOCATE( face_junction_list(cell_face)%n_external_conductors(1:3) )
    ALLOCATE( face_junction_list(cell_face)%P_matrix_list(1:3) )
    
    read(cable_model_file_unit,*)face_junction_list(cell_face)%segment_list(1:2)
    read(cable_model_file_unit,*)face_junction_list(cell_face)%n_external_conductors(1:3)
        
    n_internal=face_junction_list(cell_face)%n_internal_connection_nodes
    do face=1,3
    
      n_external=face_junction_list(cell_face)%n_external_conductors(face)
      if (n_external.ne.0) then
      
        ALLOCATE( face_junction_list(cell_face)%P_matrix_list(face)%P(1:n_internal,1:n_external) )
        do row=1,face_junction_list(cell_face)%n_internal_connection_nodes
	
          read(cable_model_file_unit,*)face_junction_list(cell_face)%P_matrix_list(face)%P(row,1:n_external)
	  
        end do
	
      end if ! n_external.gt.0
      
    end do ! next face
    
! Internal impedance stuff    
    read(cable_model_file_unit,*)face_junction_list(cell_face)%n_internal_impedance_filters
    n_filters=face_junction_list(cell_face)%n_internal_impedance_filters
    
    ALLOCATE( face_junction_list(cell_face)%Sfilter(1:n_filters) )
    ALLOCATE( face_junction_list(cell_face)%Zfilter(1:n_filters) )
    ALLOCATE( face_junction_list(cell_face)%Z_f(1:n_filters) )
   
    do filter=1,face_junction_list(cell_face)%n_internal_impedance_filters

      read(cable_model_file_unit,*)  
      CALL read_Sfilter(face_junction_list(cell_face)%Sfilter(filter),cable_model_file_unit)
      read(cable_model_file_unit,*)        
      read(cable_model_file_unit,*)face_junction_list(cell_face)%Z_f(filter)
      read(cable_model_file_unit,*)   
      CALL read_Zfilter(face_junction_list(cell_face)%Zfilter(filter),cable_model_file_unit)
      
    end do ! next filter
   
  end do ! next face junction
  
  
! READ CABLE_OUTPUT_LIST
  CALL write_line('Read cable output list',0,output_to_screen_flag)
  
  read(cable_model_file_unit,*) ! Read comment line
  
  read(cable_model_file_unit,*)n_cable_outputs
  
  ALLOCATE( cable_output_list(1:n_cable_outputs) )
  
  do output=1,n_cable_outputs
    
    read(cable_model_file_unit,*)cable_output_list(output)%cable_number
    read(cable_model_file_unit,*)cable_output_list(output)%closest_point_number
    read(cable_model_file_unit,*)cable_output_list(output)%bundle_segment_number
    read(cable_model_file_unit,*)cable_output_list(output)%output_point%cell%i, 	 &
    				     cable_output_list(output)%output_point%cell%j,	     &
    				     cable_output_list(output)%output_point%cell%k,	     &
    				     cable_output_list(output)%output_point%point
				     
    read(cable_model_file_unit,*)cable_output_list(output)%n_conductors
    
    ALLOCATE( cable_output_list(output)%conductor_list(1:cable_output_list(output)%n_conductors) )
    
    do conductor=1,cable_output_list(output)%n_conductors
      read(cable_model_file_unit,*)cable_output_list(output)%conductor_list(conductor)
    end do ! next conductor
    
  end do ! next cable output
  
  if ( .NOT.read_data_for_computation_only ) then
! read additional data for cable model checks  

  
! READ CABLE_GEOMETRY LIST
  
    read(cable_model_file_unit,*) ! Read comment line
    read(cable_model_file_unit,*)n_cable_geometries
  
    if (n_cable_geometries.NE.0) then
      ALLOCATE( cable_geometry_list(1:n_cable_geometries) )
    end if
  
    do cable_geometry=1,n_cable_geometries
  
      read(cable_model_file_unit,'(A)')cable_geometry_list(cable_geometry)%cable_geometry_type_string
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%cable_geometry_type
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_conductors
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_shielded_conductors
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_external_conductors
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_filters
      
      n_external=cable_geometry_list(cable_geometry)%n_external_conductors
      
      ALLOCATE( cable_geometry_list(cable_geometry)%external_conductor_xc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry)%external_conductor_yc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry)%external_conductor_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry)%external_dielectric_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry)%external_dielectric_permittivity(1:n_external) )

      do i=1,n_external
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_xc(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_yc(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_radius(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_radius(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_permittivity(i)
      end do
      
      n_shielded=cable_geometry_list(cable_geometry)%n_shielded_conductors
      
      ALLOCATE( cable_geometry_list(cable_geometry)%shielded_conductor_xc(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry)%shielded_conductor_yc(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry)%shielded_conductor_radius(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry)%shielded_dielectric_radius(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry)%shielded_dielectric_permittivity(1:n_shielded) )

      do i=1,n_shielded
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_conductor_xc(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_conductor_yc(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_conductor_radius(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_dielectric_radius(i)
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_dielectric_permittivity(i)
      end do
    
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%cable_offset_radius

      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_radius
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_radius
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_permittivity

! read parameters
      read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_parameters
      n_cols=cable_geometry_list(cable_geometry)%n_parameters
    
      ALLOCATE( cable_geometry_list(cable_geometry)%parameters(1:n_cols) )
    
      read(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%parameters(col),	&
                                   col=1,cable_geometry_list(cable_geometry)%n_parameters)
				   
      n_rows=cable_geometry_list(cable_geometry)%n_conductors
      n_cols=n_rows
    
      ALLOCATE( cable_geometry_list(cable_geometry)%Sc(1:n_rows) )
      ALLOCATE( cable_geometry_list(cable_geometry)%Tv(1:n_rows,1:n_cols) )
      ALLOCATE( cable_geometry_list(cable_geometry)%Ti(1:n_rows,1:n_cols) )
      ALLOCATE( cable_geometry_list(cable_geometry)%L_internal(1:n_rows,1:n_cols) )
      ALLOCATE( cable_geometry_list(cable_geometry)%C_internal(1:n_rows,1:n_cols) )
      
      ALLOCATE( cable_geometry_list(cable_geometry)%filter_number(1:n_rows,1:n_cols) )
      ALLOCATE( cable_geometry_list(cable_geometry)%Sfilter(1:cable_geometry_list(cable_geometry)%n_filters) )
      ALLOCATE( cable_geometry_list(cable_geometry)%Z_f(1:cable_geometry_list(cable_geometry)%n_filters) )
      ALLOCATE( cable_geometry_list(cable_geometry)%Zfilter(1:cable_geometry_list(cable_geometry)%n_filters) )
    
! read Sc
      read(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%Sc(row),row=1,n_rows)
    
! read Tv
      do row=1,n_rows	   
        read(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%Tv(row,col),col=1,n_cols)
      end do
    
! read Ti
      do row=1,n_rows	   
        read(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%Ti(row,col),col=1,n_cols)
      end do
    
! read L_internal
      do row=1,n_rows	   
        read(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%L_internal(row,col),col=1,n_cols)
      end do
    
! read C_internal
      do row=1,n_rows	   
        read(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%C_internal(row,col),col=1,n_cols)
      end do
 
! read filter information

! read filter_number
      do row=1,n_rows	   
        read(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%filter_number(row,col),col=1,n_cols)
      end do
    
      n_rows=cable_geometry_list(cable_geometry)%n_filters
    
! read filter information
      do row=1,n_rows
        read(cable_model_file_unit,*)     
        CALL read_Sfilter(cable_geometry_list(cable_geometry)%Sfilter(row),cable_model_file_unit)
      end do ! next row
    
      do row=1,n_rows
        read(cable_model_file_unit,*)        
        read(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%Z_f(row)
        read(cable_model_file_unit,*)         
        CALL read_Zfilter(cable_geometry_list(cable_geometry)%Zfilter(row),cable_model_file_unit)
      end do ! next row
 
    end do ! next cable geometry
  
  
! READ CABLE LIST
  
    read(cable_model_file_unit,*) ! Read comment line
    read(cable_model_file_unit,*)n_cables
  
    if (n_cables.NE.0) then
      ALLOCATE( cable_list(1:n_cables) )
    end if
    
    do cable=1,n_cables
  
      read(cable_model_file_unit,*)cable_list(cable)%cable_geometry_number
      read(cable_model_file_unit,*)cable_list(cable)%n_lines  
      
      ALLOCATE( cable_list(cable)%line_list(1:cable_list(cable)%n_lines) )
      
      read(cable_model_file_unit,*)(cable_list(cable)%line_list(i),i=1,cable_list(cable)%n_lines) 
      read(cable_model_file_unit,*)cable_list(cable)%junction_1  
      read(cable_model_file_unit,*)cable_list(cable)%junction_2  
      read(cable_model_file_unit,*)cable_list(cable)%number_of_cable_segments  
      
      ALLOCATE( cable_list(cable)%cable_segment_list(1:cable_list(cable)%number_of_cable_segments) )
      ALLOCATE( cable_list(cable)%direction_sign_list(1:cable_list(cable)%number_of_cable_segments) )
      ALLOCATE( cable_list(cable)%bundle_segment_list(1:cable_list(cable)%number_of_cable_segments) )
      
      read(cable_model_file_unit,*)(cable_list(cable)%cable_segment_list(i),i=1,cable_list(cable)%number_of_cable_segments)
      read(cable_model_file_unit,*)(cable_list(cable)%direction_sign_list(i),i=1,cable_list(cable)%number_of_cable_segments)
      read(cable_model_file_unit,*)(cable_list(cable)%bundle_segment_list(i),i=1,cable_list(cable)%number_of_cable_segments)  
  
    end do ! next cable

  end if ! .NOT. read_data_for_computation_only
  
  CALL close_file(cable_model_file_unit)

  CALL write_line('FINISHED: read_cable_model',0,output_to_screen_flag)
    
  RETURN
  
9000 CALL write_line('ERROR in read_cable model',0,.TRUE.)
     CALL write_line('Number of segments in a cell centre junction should be 6',0,.TRUE.)
     CALL write_line_integer('cell_centre_junction_list(cell)%n_segments',	&
                              cell_centre_junction_list(cell)%n_segments,0,.TRUE.)
     STOP
  
9010 CALL write_line('ERROR in read_cable model',0,.TRUE.)
     CALL write_line('Number of segments in a face junction should be 2',0,.TRUE.)
     CALL write_line_integer('face_junction_list(cell_face)%n_segments',	&
                              face_junction_list(cell_face)%n_segments,0,.TRUE.)
     CALL write_line_integer('face_junction number',cell_face,0,.TRUE.)
     STOP
  
  
END SUBROUTINE read_cable_model
