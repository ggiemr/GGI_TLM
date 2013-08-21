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
!SUBROUTINE write_cable_model
!
! NAME
!     SUBROUTINE write_cable_model
!
! DESCRIPTION
!      write the cable model data to file
!      The first set of data written is that required for computation in GGI_TLM
!      Following this we write additional data which is read by GGI_TLM_cable_model_checks
!      in order to examine the detail of the construction of the cable model
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2012 CJS
!
!
SUBROUTINE write_cable_model()

USE TLM_general
USE Cables
USE File_information
USE filter_types
USE filter_functions

IMPLICIT NONE

! local variables

  integer 	:: segment
  integer 	:: segment_geometry
  integer	:: cell
  integer	:: cable
  integer	:: cable_geometry
  integer	:: output
  integer	:: conductor
  integer	:: n_external
  integer	:: face
  integer	:: row,col
  integer	:: n_rows,n_cols
  integer	:: i
  integer 	:: filter
  
! START

  CALL write_line('CALLED: write_cable_model',0,output_to_screen_flag)

! Open cable model file

  CALL open_file(cable_model_file_unit,cable_model_file_extn)
  
! WRITE BUNDLE_SEGMENT_LIST
  write(cable_model_file_unit,*)'BUNDLE SEGMENT LIST'
  write(cable_model_file_unit,*)n_bundle_segments,' ! NUMBER OF BUNDLE SEGMENTS'
  do segment=1,n_bundle_segments
    
    write(cable_model_file_unit,8000)bundle_segment_list(segment)%cable_segment%segment_point(1)%cell%i,	&
    				     bundle_segment_list(segment)%cable_segment%segment_point(1)%cell%j,	&
    				     bundle_segment_list(segment)%cable_segment%segment_point(1)%cell%k,	&
    				     bundle_segment_list(segment)%cable_segment%segment_point(1)%point,	&
    				     bundle_segment_list(segment)%cable_segment%segment_point(2)%cell%i,	&
    				     bundle_segment_list(segment)%cable_segment%segment_point(2)%cell%j,	&
    				     bundle_segment_list(segment)%cable_segment%segment_point(2)%cell%k,	&
    				     bundle_segment_list(segment)%cable_segment%segment_point(2)%point,	&
				     ' ! cell_point1, cell_point2'
    
    write(cable_model_file_unit,*)bundle_segment_list(segment)%n_cables,' ! n_cables'
    
    do cable=1,bundle_segment_list(segment)%n_cables    
      write(cable_model_file_unit,*)bundle_segment_list(segment)%cable_list(cable),' ! cable number'
    end do ! next cable
    
    write(cable_model_file_unit,*)bundle_segment_list(segment)%n_conductors,' ! n_conductors'
    
    write(cable_model_file_unit,*)bundle_segment_list(segment)%n_filters,' ! n_filters'
    
    do row=1,bundle_segment_list(segment)%n_conductors
      write(cable_model_file_unit,*)bundle_segment_list(segment)%direction_sign_list(row),' ! direction sign'
    end do ! next conductor

    write(cable_model_file_unit,*)bundle_segment_list(segment)%bundle_segment_geometry,' ! bundle_segment_geometry'
    
    do row=1,bundle_segment_list(segment)%n_conductors
      do col=1,bundle_segment_list(segment)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_list(segment)%L(row,col),' ! L(i,j)'          
      end do ! next col
    end do ! next row
     
    do row=1,bundle_segment_list(segment)%n_conductors
      do col=1,bundle_segment_list(segment)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_list(segment)%C(row,col),' ! C(i,j)'       
      end do ! next col
    end do ! next row
   
    do row=1,bundle_segment_list(segment)%n_conductors
      do col=1,bundle_segment_list(segment)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_list(segment)%R(row,col),' ! R(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_list(segment)%n_conductors
      do col=1,bundle_segment_list(segment)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_list(segment)%Tv(row,col),' ! Tv(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_list(segment)%n_conductors
      do col=1,bundle_segment_list(segment)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_list(segment)%Ti(row,col),' ! Ti(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_list(segment)%n_conductors
      write(cable_model_file_unit,*)bundle_segment_list(segment)%SC(row),' ! SC(i)'         
    end do ! next row
    
    do row=1,bundle_segment_list(segment)%n_conductors
      write(cable_model_file_unit,*)bundle_segment_list(segment)%excitation_function(row),' ! excitation_function(i)'          
    end do ! next row
    
    do row=1,bundle_segment_list(segment)%n_conductors
      do col=1,bundle_segment_list(segment)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_list(segment)%filter_number(row,col),' ! filter_number(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_list(segment)%n_filters
      write(cable_model_file_unit,*)'! Sfilter number',row          
      CALL write_Sfilter(bundle_segment_list(segment)%Sfilter(row),cable_model_file_unit)
    end do ! next row
    
    do row=1,bundle_segment_list(segment)%n_filters
      write(cable_model_file_unit,*)'! Z_f',row          
      write(cable_model_file_unit,*)bundle_segment_list(segment)%Z_f(row)
      write(cable_model_file_unit,*)'! Zfilter number',row          
      CALL write_Zfilter(bundle_segment_list(segment)%Zfilter(row),cable_model_file_unit)
    end do ! next row
 	
  end do ! next bundle segment
  
! WRITE BUNDLE_SEGMENT_GEOMETYRY LIST
  write(cable_model_file_unit,*)'BUNDLE SEGMENT GEOMETRY LIST'
  write(cable_model_file_unit,*)n_bundle_segment_geometries,' ! NUMBER OF BUNDLE SEGMENT GEOMETRIES'
  do segment_geometry=1,n_bundle_segment_geometries
        
    write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%n_cables,' ! n_cables'
    
    do cable=1,bundle_segment_geometry_list(segment_geometry)%n_cables    
      write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%cable_list(cable),' ! cable number'
    end do ! next cable
    
    write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%n_conductors,' ! n_conductors'
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%xc(row),	&
                                    bundle_segment_geometry_list(segment_geometry)%yc(row),	&
                                    bundle_segment_geometry_list(segment_geometry)%rc(row),	&
                                    bundle_segment_geometry_list(segment_geometry)%ri(row),' conductor x, y, r, ri'       
    end do ! next row   
    
    write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%cable_bundle_radius,&
                                 ' ! cable_bundle_radius'
    write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%TLM_cell_equivalent_radius,&
                                 ' ! TLM_cell_equivalent_radius'
    write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%TLM_reference_radius_rL,&
                                 ' ! TLM_reference_radius_rL'
    write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%TLM_reference_radius_rC,&
                                 ' ! TLM_reference_radius_rC'
      
    write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%n_filters,' ! n_filters'
   
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%L(row,col),' ! L(i,j)'          
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%C(row,col),' ! C(i,j)'       
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%R(row,col),' ! R(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Zlink(row,col),' ! Zlink(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Ylink(row,col),' ! Ylink(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%ZLstub(row,col),' ! ZLstub(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Yf(row,col),' ! Yf(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Tv(row,col),' ! Tv(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Ti(row,col),' ! Ti(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%SC(row),' ! SC(i)'         
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
      do col=1,bundle_segment_geometry_list(segment_geometry)%n_conductors
        write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%filter_number(row,col),' ! filter_number(i,j)'         
      end do ! next col
    end do ! next row
    
    do row=1,bundle_segment_geometry_list(segment_geometry)%n_filters
      write(cable_model_file_unit,*)'! Sfilter number',row          
      CALL write_Sfilter(bundle_segment_geometry_list(segment_geometry)%Sfilter(row),cable_model_file_unit)
    end do ! next row
    
    do row=1,bundle_segment_list(segment_geometry)%n_filters
      write(cable_model_file_unit,*)'! Z_f',row          
      write(cable_model_file_unit,*)bundle_segment_geometry_list(segment_geometry)%Z_f(row)
      write(cable_model_file_unit,*)'! Zfilter number',row          
      CALL write_Zfilter(bundle_segment_geometry_list(segment_geometry)%Zfilter(row),cable_model_file_unit)
    end do ! next row
	
  end do ! next bundle segment_geometry

! WRITE CELL_CENTRE_JUNCTION_LIST
  write(cable_model_file_unit,*)'CELL_CENTRE_JUNCTION_LIST'
  write(cable_model_file_unit,*)n_cell_centre_junctions,' ! NUMBER OF CELL CENTRE JUNCTIONS'

  do cell=1,n_cell_centre_junctions
  
    write(cable_model_file_unit,8010)cell_centre_junction_list(cell)%cell_point%cell%i,	&
                                     cell_centre_junction_list(cell)%cell_point%cell%j,	&
                                     cell_centre_junction_list(cell)%cell_point%cell%k,	&
                                     cell_centre_junction_list(cell)%cell_point%point,' ! ix,iy,iz,face'
				     
    write(cable_model_file_unit,8020)cell_centre_junction_list(cell)%n_internal_connection_nodes,' ! n_internal_connection_nodes'
    write(cable_model_file_unit,8020)cell_centre_junction_list(cell)%n_segments,' ! n_segments '
    write(cable_model_file_unit,8030)cell_centre_junction_list(cell)%segment_list(1:6),' ! segment data (1:6)'
    write(cable_model_file_unit,8035)cell_centre_junction_list(cell)%n_external_conductors(1:7),' ! n_external conductor data (1:7)'

    write(cable_model_file_unit,8020)cell_centre_junction_list(cell)%number_of_cable_junctions,' !number_of_cable_junctions '

    do i=1,cell_centre_junction_list(cell)%number_of_cable_junctions
      write(cable_model_file_unit,8020)cell_centre_junction_list(cell)%cable_junction_list(i)
    end do
    
    do face=1,7
    
      n_external=cell_centre_junction_list(cell)%n_external_conductors(face)
      if (n_external.ne.0) then
      
        do row=1,cell_centre_junction_list(cell)%n_internal_connection_nodes
	
          write(cable_model_file_unit,8040)cell_centre_junction_list(cell)%P_matrix_list(face)%P(row,1:n_external)
	  
        end do
	
      end if ! n_external.gt.0
      
    end do ! next face
    
! Internal impedance stuff    
    write(cable_model_file_unit,8020)cell_centre_junction_list(cell)%n_internal_impedance_filters,	&
                                     ' !number_of_internal impedance filters '
   
    do filter=1,cell_centre_junction_list(cell)%n_internal_impedance_filters

      write(cable_model_file_unit,*)'! Sfilter number',filter         
      CALL write_Sfilter(cell_centre_junction_list(cell)%Sfilter(filter),cable_model_file_unit)
      write(cable_model_file_unit,*)'! Z_f',filter          
      write(cable_model_file_unit,*)cell_centre_junction_list(cell)%Z_f(filter)
      write(cable_model_file_unit,*)'! Zfilter number',filter      
      CALL write_Zfilter(cell_centre_junction_list(cell)%Zfilter(filter),cable_model_file_unit)
      
    end do ! next filter
   
  end do ! next cable cell

! WRITE FACE_JUNCTION_LIST
  write(cable_model_file_unit,*)'FACE_JUNCTION_LIST'
  write(cable_model_file_unit,*)n_face_junctions,' ! NUMBER OF FACE JUNCTIONS'

  do cell=1,n_face_junctions
  
    write(cable_model_file_unit,8010)face_junction_list(cell)%cell_point%cell%i,	&
                                     face_junction_list(cell)%cell_point%cell%j,	&
                                     face_junction_list(cell)%cell_point%cell%k,	&
                                     face_junction_list(cell)%cell_point%point,' ! ix,iy,iz,face'
				     
    write(cable_model_file_unit,8020)face_junction_list(cell)%n_internal_connection_nodes,' ! n_internal_connection_nodes'
    do row=1,face_junction_list(cell)%n_internal_connection_nodes
      write(cable_model_file_unit,8020)face_junction_list(cell)%BC(row),' ! Boundary condition data '
    end do
    write(cable_model_file_unit,8020)face_junction_list(cell)%n_segments,' ! n_segments '
    write(cable_model_file_unit,8050)face_junction_list(cell)%segment_list(1:2),' ! segment data (1:2)'
    write(cable_model_file_unit,8055)face_junction_list(cell)%n_external_conductors(1:3),' ! n_external conductor data (1:3)'
    
    do face=1,3
    
      n_external=face_junction_list(cell)%n_external_conductors(face)
      if (n_external.ne.0) then
      
        do row=1,face_junction_list(cell)%n_internal_connection_nodes
	
          write(cable_model_file_unit,8040)face_junction_list(cell)%P_matrix_list(face)%P(row,1:n_external)
	  
        end do
	
      end if ! n_external.gt.0
      
    end do ! next face
    
! Internal impedance stuff    
    write(cable_model_file_unit,8020)face_junction_list(cell)%n_internal_impedance_filters,	&
                                     ' !number_of_internal impedance filters '
   
    do filter=1,face_junction_list(cell)%n_internal_impedance_filters

      write(cable_model_file_unit,*)'! Sfilter number',filter         
      CALL write_Sfilter(face_junction_list(cell)%Sfilter(filter),cable_model_file_unit)
      write(cable_model_file_unit,*)'! Z_f',filter          
      write(cable_model_file_unit,*)face_junction_list(cell)%Z_f(filter)
      write(cable_model_file_unit,*)'! Zfilter number',filter      
      CALL write_Zfilter(face_junction_list(cell)%Zfilter(filter),cable_model_file_unit)
      
    end do ! next filter
   
  end do ! next face junction

! WRITE CABLE_OUTPUT_LIST
  write(cable_model_file_unit,*)'CABLE_OUTPUT_LIST'
  write(cable_model_file_unit,*)n_cable_outputs,' ! NUMBER OF CABLE OUTPUTS'
  
  do output=1,n_cable_outputs
    
    write(cable_model_file_unit,*)cable_output_list(output)%cable_number,         ' ! cable_number'
    write(cable_model_file_unit,*)cable_output_list(output)%closest_point_number, ' !  closest_point_number'
    write(cable_model_file_unit,*)cable_output_list(output)%bundle_segment_number,' ! bundle_segment_number '
    write(cable_model_file_unit,8060)cable_output_list(output)%output_point%cell%i,	     &
    				     cable_output_list(output)%output_point%cell%j,	     &
    				     cable_output_list(output)%output_point%cell%k,	     &
    				     cable_output_list(output)%output_point%point, &
				     ' ! cell_point'
    write(cable_model_file_unit,*)cable_output_list(output)%n_conductors,' ! n_conductors '
    
    do conductor=1,cable_output_list(output)%n_conductors
      write(cable_model_file_unit,*)cable_output_list(output)%conductor_list(conductor),' ! conductor list '
    end do ! next conductor
    
  end do ! next cable output
  
! WRITE CABLE_GEOMETRY LIST
  write(cable_model_file_unit,*)'CABLE_GEOMETRY_LIST'
  write(cable_model_file_unit,*)n_cable_geometries,' ! NUMBER OF CABLE GEOMETRIES'
  
  do cable_geometry=1,n_cable_geometries
    
    write(cable_model_file_unit,*)trim(cable_geometry_list(cable_geometry)%cable_geometry_type_string)
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%cable_geometry_type
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_conductors
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_shielded_conductors
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_external_conductors
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_filters

    do i=1,cable_geometry_list(cable_geometry)%n_external_conductors
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_xc(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_yc(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_radius(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_radius(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_permittivity(i)
    end do

    do i=1,cable_geometry_list(cable_geometry)%n_shielded_conductors
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_conductor_xc(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_conductor_yc(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_conductor_radius(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_dielectric_radius(i)
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%shielded_dielectric_permittivity(i)
    end do
    
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%cable_offset_radius

    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_conductor_radius
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_radius
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%external_dielectric_permittivity

! write parameters
    write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%n_parameters
    
    n_cols=cable_geometry_list(cable_geometry)%n_parameters
    write(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%parameters(col),	&
                                   col=1,cable_geometry_list(cable_geometry)%n_parameters)
				   
    n_rows=cable_geometry_list(cable_geometry)%n_conductors
    n_cols=n_rows
    
! write Sc
    write(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%Sc(row),row=1,n_rows)
    
! write Tv
    do row=1,n_rows	   
      write(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%Tv(row,col),col=1,n_cols)
    end do
    
! write Ti
    do row=1,n_rows	   
      write(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%Ti(row,col),col=1,n_cols)
    end do
    
! write L_internal
    do row=1,n_rows	   
      write(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%L_internal(row,col),col=1,n_cols)
    end do
    
! write C_internal
    do row=1,n_rows	   
      write(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%C_internal(row,col),col=1,n_cols)
    end do
    
! write filter information

! write filter_number
    do row=1,n_rows	   
      write(cable_model_file_unit,*)(cable_geometry_list(cable_geometry)%filter_number(row,col),col=1,n_cols)
    end do
    
    n_rows=cable_geometry_list(cable_geometry)%n_filters
    
! write filter information
    do row=1,n_rows
      write(cable_model_file_unit,*)'! Sfilter number',row          
      CALL write_Sfilter(cable_geometry_list(cable_geometry)%Sfilter(row),cable_model_file_unit)
    end do ! next row
    
    do row=1,n_rows
      write(cable_model_file_unit,*)'! Z_f',row          
      write(cable_model_file_unit,*)cable_geometry_list(cable_geometry)%Z_f(row)
      write(cable_model_file_unit,*)'! Zfilter number',row          
      CALL write_Zfilter(cable_geometry_list(cable_geometry)%Zfilter(row),cable_model_file_unit)
    end do ! next row
 
  end do ! next cable geometry
    
! WRITE CABLE LIST
  write(cable_model_file_unit,*)'CABLE_LIST'
  write(cable_model_file_unit,*)n_cables,' ! NUMBER OF CABLES'
  do cable=1,n_cables
  
    write(cable_model_file_unit,*)cable_list(cable)%cable_geometry_number
    write(cable_model_file_unit,*)cable_list(cable)%n_lines  
    write(cable_model_file_unit,*)(cable_list(cable)%line_list(i),i=1,cable_list(cable)%n_lines) 
    write(cable_model_file_unit,*)cable_list(cable)%junction_1  
    write(cable_model_file_unit,*)cable_list(cable)%junction_2  
    write(cable_model_file_unit,*)cable_list(cable)%number_of_cable_segments  
    write(cable_model_file_unit,*)(cable_list(cable)%cable_segment_list(i),	&
                                   i=1,cable_list(cable)%number_of_cable_segments)
    write(cable_model_file_unit,*)(cable_list(cable)%direction_sign_list(i),	&
                                   i=1,cable_list(cable)%number_of_cable_segments)
    write(cable_model_file_unit,*)(cable_list(cable)%bundle_segment_list(i),	&
                                   i=1,cable_list(cable)%number_of_cable_segments)  
  
  end do ! next cable
  
  CALL close_file(cable_model_file_unit)

  CALL write_line('FINISHED: write_cable_model',0,output_to_screen_flag)
    
  RETURN
     
8000 format(8I8,A)    
8010 format(4I8,A) 
8020 format(I8,A) 
8030 format(6I8,A) 
8035 format(7I8,A) 
8040 format(1000I3) 
8050 format(2I8,A) 
8055 format(3I8,A) 
8060 format(4I8,A)    
  
END SUBROUTINE write_cable_model
