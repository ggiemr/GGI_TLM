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
! SUBROUTINE deallocate_geometry
! SUBROUTINE deallocate_materials
! SUBROUTINE deallocate_cables
! SUBROUTINE deallocate_outputs
! SUBROUTINE deallocate_excitations
! SUBROUTINE deallocate_mode_stir
! SUBROUTINE deallocate_mesh
!
! NAME
!     SUBROUTINE deallocate_geometry
!
! DESCRIPTION
!     deallocate_geometry:
!
!     deallocate all allocatable structures related to geometry definition
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE deallocate_geometry

USE TLM_general
USE geometry

IMPLICIT NONE

! local variables

integer	:: volume_number
integer	:: surface_number

! START

  
  CALL write_line('CALLED: deallocate_geometry',0,output_to_screen_flag)
  
  if ( allocated( problem_volumes ) ) then
  
    do volume_number=1,n_volumes
  
      if ( allocated( problem_volumes(volume_number)%tet_list ) ) then
        DEALLOCATE( problem_volumes(volume_number)%tet_list )
      end if
  
      if ( allocated( problem_volumes(volume_number)%cell_list ) ) then
        DEALLOCATE( problem_volumes(volume_number)%cell_list )
      end if
    
    end do ! next volume number
    
    DEALLOCATE( problem_volumes )
    
  end if
  
  if ( allocated( problem_surfaces ) ) then
  
    do surface_number=1,n_surfaces
  
      if ( allocated( problem_surfaces(surface_number)%triangle_list ) ) then
        DEALLOCATE( problem_surfaces(surface_number)%triangle_list )
      end if
  
      if ( allocated( problem_surfaces(surface_number)%face_list ) ) then
        DEALLOCATE( problem_surfaces(surface_number)%face_list )
      end if
    
    end do ! next surface number
    
    DEALLOCATE( problem_surfaces )
    
  end if
  
  if ( allocated( problem_points ) ) then
      
    DEALLOCATE( problem_points )
    
  end if
  
  CALL write_line('FINISHED: deallocate_geometry',0,output_to_screen_flag)

  RETURN
  
END SUBROUTINE deallocate_geometry
! 
!
! NAME
!     SUBROUTINE deallocate_materials
!
! DESCRIPTION
!     deallocate_materials:
!
!     deallocate all allocatable structures related to geometry definition
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE deallocate_materials

USE TLM_general
USE TLM_volume_materials
USE TLM_surface_materials

IMPLICIT NONE

! local variables

integer	:: i,j

! START

  CALL write_line('CALLED: deallocate_materials',0,output_to_screen_flag)
  
  if (allocated( volume_material_list )) then 
    do i=1,n_volume_materials
    
      CALL deallocate_Sfilter( volume_material_list(i)%eps_S )
      CALL deallocate_Sfilter( volume_material_list(i)%Zs_eps_S )
      CALL deallocate_Zfilter( volume_material_list(i)%Zs_eps_Z )
      CALL deallocate_Sfilter( volume_material_list(i)%mu_S )
      CALL deallocate_Sfilter( volume_material_list(i)%Zs_mu_S )
      CALL deallocate_Zfilter( volume_material_list(i)%Zs_mu_Z )

      if ( allocated(volume_material_list(i)%volume_list) ) then
        DEALLOCATE( volume_material_list(i)%volume_list )
      end if
      
    end do ! next material in list
  end if 
  
  if (allocated( volume_material_Zs_eps_filter_data )) then 
    do i=1,volume_material_storage
      CALL deallocate_Zfilter_data( volume_material_Zs_eps_filter_data(i) )
    end do 
    DEALLOCATE( volume_material_Zs_eps_filter_data )
  end if 
  
  if (allocated( volume_material_Zs_mu_filter_data )) then 
    do i=1,volume_material_storage
      CALL deallocate_Zfilter_data( volume_material_Zs_mu_filter_data(i) )     
    end do
    DEALLOCATE( volume_material_Zs_mu_filter_data )
  end if 
  
  
  if (allocated( surface_material_list )) then 
    do i=1,n_surface_materials
    
      CALL deallocate_Sfilter( surface_material_list(i)%Z11_S )
      CALL deallocate_Zfilter( surface_material_list(i)%Z11_Z )
      CALL deallocate_Sfilter( surface_material_list(i)%Z12_S )
      CALL deallocate_Zfilter( surface_material_list(i)%Z12_Z )
      CALL deallocate_Sfilter( surface_material_list(i)%Z21_S )
      CALL deallocate_Zfilter( surface_material_list(i)%Z21_Z )
      CALL deallocate_Sfilter( surface_material_list(i)%Z22_S )
      CALL deallocate_Zfilter( surface_material_list(i)%Z22_Z )
      
      if ( allocated(surface_material_list(i)%surface_list) ) then
	DEALLOCATE( surface_material_list(i)%surface_list )
	DEALLOCATE( surface_material_list(i)%surface_orientation_list )
      end if
      
    end do ! next material in list
  end if 
  
  if (allocated( surface_material_Z11_filter_data )) then 
    do i=1,surface_material_storage
      CALL deallocate_Zfilter_data( surface_material_Z11_filter_data(i) )
    end do 
    DEALLOCATE( surface_material_Z11_filter_data )
  end if 
  
  if (allocated( surface_material_Z12_filter_data )) then 
    do i=1,surface_material_storage
      CALL deallocate_Zfilter_data( surface_material_Z12_filter_data(i) )     
    end do 
    DEALLOCATE( surface_material_Z12_filter_data )
  end if 
  
  if (allocated( surface_material_Z21_filter_data )) then 
    do i=1,surface_material_storage
      CALL deallocate_Zfilter_data( surface_material_Z21_filter_data(i) )
    end do 
    DEALLOCATE( surface_material_Z21_filter_data )
  end if 
  
  if (allocated( surface_material_Z22_filter_data )) then 
    do i=1,surface_material_storage
      CALL deallocate_Zfilter_data( surface_material_Z22_filter_data(i) )     
    end do 
    DEALLOCATE( surface_material_Z22_filter_data )
  end if 
  
  CALL write_line('FINISHED: deallocate_materials',0,output_to_screen_flag)

  RETURN
  
END SUBROUTINE deallocate_materials
! 
!
! NAME
!     SUBROUTINE deallocate_cables
!
! DESCRIPTION
!     deallocate_cables:
!
!     deallocate all allocatable structures related to geometry definition
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE deallocate_cables

USE TLM_general
USE Cables

IMPLICIT NONE

! local variables

integer	:: i,p

! START

  CALL write_line('CALLED: deallocate_cables',0,output_to_screen_flag)
  
  if (allocated( cable_geometry_list )) then
  
    do i=1,n_cable_geometries
    
      if(allocated( cable_geometry_list(i)%external_conductor_xc ))  		&
        DEALLOCATE( cable_geometry_list(i)%external_conductor_xc ) 
      if(allocated( cable_geometry_list(i)%external_conductor_yc ))		&
        DEALLOCATE( cable_geometry_list(i)%external_conductor_yc ) 
      if(allocated( cable_geometry_list(i)%external_conductor_radius ))		&
        DEALLOCATE( cable_geometry_list(i)%external_conductor_radius ) 
      if(allocated( cable_geometry_list(i)%external_dielectric_radius ))	&
        DEALLOCATE( cable_geometry_list(i)%external_dielectric_radius ) 
      if(allocated( cable_geometry_list(i)%external_dielectric_permittivity ))	&
        DEALLOCATE( cable_geometry_list(i)%external_dielectric_permittivity ) 
    
      if(allocated( cable_geometry_list(i)%shielded_conductor_xc ))  		&
        DEALLOCATE( cable_geometry_list(i)%shielded_conductor_xc ) 
      if(allocated( cable_geometry_list(i)%shielded_conductor_yc ))		&
        DEALLOCATE( cable_geometry_list(i)%shielded_conductor_yc ) 
      if(allocated( cable_geometry_list(i)%shielded_conductor_radius ))		&
        DEALLOCATE( cable_geometry_list(i)%shielded_conductor_radius ) 
      if(allocated( cable_geometry_list(i)%shielded_dielectric_radius ))	&
        DEALLOCATE( cable_geometry_list(i)%shielded_dielectric_radius ) 
      if(allocated( cable_geometry_list(i)%shielded_dielectric_permittivity ))	&
        DEALLOCATE( cable_geometry_list(i)%shielded_dielectric_permittivity ) 

      if(allocated( cable_geometry_list(i)%L_internal ))  DEALLOCATE( cable_geometry_list(i)%L_internal ) 
      if(allocated( cable_geometry_list(i)%C_internal ))  DEALLOCATE( cable_geometry_list(i)%C_internal ) 
      if(allocated( cable_geometry_list(i)%R_internal ))  DEALLOCATE( cable_geometry_list(i)%R_internal ) 
      if(allocated( cable_geometry_list(i)%Tv )) DEALLOCATE( cable_geometry_list(i)%Tv ) 
      if(allocated( cable_geometry_list(i)%Ti )) DEALLOCATE( cable_geometry_list(i)%Ti ) 
      if(allocated( cable_geometry_list(i)%SC )) DEALLOCATE( cable_geometry_list(i)%SC ) 
      if(allocated( cable_geometry_list(i)%parameters ))  DEALLOCATE( cable_geometry_list(i)%parameters ) 

      if(allocated( cable_geometry_list(i)%filter_number ))  DEALLOCATE( cable_geometry_list(i)%filter_number ) 
      if(allocated( cable_geometry_list(i)%Sfilter ))  DEALLOCATE( cable_geometry_list(i)%Sfilter ) 
      if(allocated( cable_geometry_list(i)%Zfilter ))  DEALLOCATE( cable_geometry_list(i)%Zfilter ) 
      if(allocated( cable_geometry_list(i)%Z_f ))  DEALLOCATE( cable_geometry_list(i)%Z_f ) 

    end do
  
    DEALLOCATE( cable_geometry_list )
  
  end if ! cable_geometry_list
  
  
  if (allocated( cable_list )) then
  
    do i=1,n_cables
      if(allocated( cable_list(i)%line_list ))           DEALLOCATE( cable_list(i)%line_list ) 
      if(allocated( cable_list(i)%cable_segment_list ))  DEALLOCATE( cable_list(i)%cable_segment_list ) 
      if(allocated( cable_list(i)%direction_sign_list ))  DEALLOCATE( cable_list(i)%direction_sign_list ) 
      if(allocated( cable_list(i)%bundle_segment_list )) DEALLOCATE( cable_list(i)%bundle_segment_list ) 
    end do
  
    DEALLOCATE( cable_list )
  
  end if ! cable_list
  
  
  if (allocated( cable_junction_list )) then
  
    do i=1,n_cable_junctions
      if(allocated( cable_junction_list(i)%cable_list )) &
        DEALLOCATE( cable_junction_list(i)%cable_list ) 
      if(allocated( cable_junction_list(i)%cable_end_list )) &
        DEALLOCATE( cable_junction_list(i)%cable_end_list  ) 
      if(allocated( cable_junction_list(i)%n_external_conductors )) &
        DEALLOCATE( cable_junction_list(i)%n_external_conductors ) 
      if(allocated( cable_junction_list(i)%excitation_function )) &
        DEALLOCATE( cable_junction_list(i)%excitation_function  ) 
      if(allocated( cable_junction_list(i)%resistance )) &
        DEALLOCATE( cable_junction_list(i)%resistance ) 
	
      if(allocated( cable_junction_list(i)%Pmatrix )) then
        do p=1,cable_junction_list(i)%number_of_cables
          if(allocated( cable_junction_list(i)%Pmatrix(p)%P )) &
	    DEALLOCATE( cable_junction_list(i)%Pmatrix(p)%P )
        end do
        DEALLOCATE( cable_junction_list(i)%Pmatrix ) 
      end if ! Pmatrix is allocated
      
      if(allocated( cable_junction_list(i)%BC )) &
        DEALLOCATE( cable_junction_list(i)%BC ) 
      if(allocated( cable_junction_list(i)%node_1 )) &
        DEALLOCATE( cable_junction_list(i)%node_1 ) 
      if(allocated( cable_junction_list(i)%node_2 )) &
        DEALLOCATE( cable_junction_list(i)%node_2 ) 
      if(allocated( cable_junction_list(i)%Sfilter )) &
        DEALLOCATE( cable_junction_list(i)%Sfilter ) 
  
    end do
  
    DEALLOCATE( cable_junction_list )
  
  end if ! cable_junction_list
  
  
  if (allocated( cell_centre_junction_list )) then
  
    do i=1,n_cell_centre_junctions
    
      if(allocated( cell_centre_junction_list(i)%segment_list )) &
        DEALLOCATE( cell_centre_junction_list(i)%segment_list ) 
      if(allocated( cell_centre_junction_list(i)%n_external_conductors )) &
        DEALLOCATE( cell_centre_junction_list(i)%n_external_conductors ) 
      if(allocated( cell_centre_junction_list(i)%BC )) &
        DEALLOCATE( cell_centre_junction_list(i)%BC ) 
      if(allocated( cell_centre_junction_list(i)%P_matrix_list )) then
        do p=1,7
          if(allocated( cell_centre_junction_list(i)%P_matrix_list(p)%P )) &
	    DEALLOCATE( cell_centre_junction_list(i)%P_matrix_list(p)%P )
        end do
        DEALLOCATE( cell_centre_junction_list(i)%P_matrix_list ) 
      end if ! Pmatrix is allocated
      
      if( allocated( cell_centre_junction_list(i)%Sfilter ) ) 	&
        DEALLOCATE( cell_centre_junction_list(i)%Sfilter )
      if( allocated( cell_centre_junction_list(i)%Zfilter ) ) 	&
        DEALLOCATE( cell_centre_junction_list(i)%Zfilter )      
      if( allocated( cell_centre_junction_list(i)%Z_f ) ) 	&
        DEALLOCATE( cell_centre_junction_list(i)%Z_f )      
      if( allocated( cell_centre_junction_list(i)%Yf ) ) 	&
        DEALLOCATE( cell_centre_junction_list(i)%Yf )
    
    end do ! next cell centre junction
  
    DEALLOCATE( cell_centre_junction_list )
  
  end if ! cell_centre_junction_list
  
  
  if (allocated( face_junction_list )) then
  
    do i=1,n_face_junctions
    
      if(allocated( face_junction_list(i)%segment_list )) &
        DEALLOCATE( face_junction_list(i)%segment_list ) 
      if(allocated( face_junction_list(i)%n_external_conductors )) &
        DEALLOCATE( face_junction_list(i)%n_external_conductors ) 
      if(allocated( face_junction_list(i)%BC )) &
        DEALLOCATE( face_junction_list(i)%BC ) 
      if(allocated( face_junction_list(i)%P_matrix_list )) then
        do p=1,face_junction_list(i)%n_segments
          if(allocated( face_junction_list(i)%P_matrix_list(p)%P )) &
	    DEALLOCATE( face_junction_list(i)%P_matrix_list(p)%P )
        end do
        DEALLOCATE( face_junction_list(i)%P_matrix_list ) 
      end if ! Pmatrix is allocated
    
    end do ! next face junction
  
    DEALLOCATE( face_junction_list )
  
  end if ! face_junction_list
  
  if (allocated( bundle_segment_list )) then
  
    do i=1,n_bundle_segments
      if(allocated( bundle_segment_list(i)%cable_list )) &
        DEALLOCATE( bundle_segment_list(i)%cable_list ) 
      if(allocated( bundle_segment_list(i)%direction_sign_list )) &
        DEALLOCATE( bundle_segment_list(i)%direction_sign_list ) 
      if(allocated( bundle_segment_list(i)%L )) &
        DEALLOCATE( bundle_segment_list(i)%L ) 
      if(allocated( bundle_segment_list(i)%C )) &
        DEALLOCATE( bundle_segment_list(i)%C ) 
      if(allocated( bundle_segment_list(i)%R )) &
        DEALLOCATE( bundle_segment_list(i)%R ) 
      if(allocated( bundle_segment_list(i)%Tv )) &
        DEALLOCATE( bundle_segment_list(i)%Tv ) 
      if(allocated( bundle_segment_list(i)%Ti )) &
        DEALLOCATE( bundle_segment_list(i)%Ti ) 
      if(allocated( bundle_segment_list(i)%SC )) &
        DEALLOCATE( bundle_segment_list(i)%SC ) 
      if(allocated( bundle_segment_list(i)%Zlink )) &
        DEALLOCATE( bundle_segment_list(i)%Zlink ) 
      if(allocated( bundle_segment_list(i)%Ylink )) &
        DEALLOCATE( bundle_segment_list(i)%Ylink ) 
      if(allocated( bundle_segment_list(i)%ZLstub )) &
        DEALLOCATE( bundle_segment_list(i)%ZLstub ) 
      if(allocated( bundle_segment_list(i)%Yf )) &
        DEALLOCATE( bundle_segment_list(i)%Yf ) 
      if(allocated( bundle_segment_list(i)%Vlink )) &
        DEALLOCATE( bundle_segment_list(i)%Vlink ) 
      if(allocated( bundle_segment_list(i)%VLstub )) &
        DEALLOCATE( bundle_segment_list(i)%VLstub ) 
      if(allocated( bundle_segment_list(i)%Vsource )) &
        DEALLOCATE( bundle_segment_list(i)%Vsource ) 
      if(allocated( bundle_segment_list(i)%Iw_centre )) &
        DEALLOCATE( bundle_segment_list(i)%Iw_centre ) 
      if(allocated( bundle_segment_list(i)%Iw_face )) &
        DEALLOCATE( bundle_segment_list(i)%Iw_face ) 
      if(allocated( bundle_segment_list(i)%excitation_function )) &
        DEALLOCATE( bundle_segment_list(i)%excitation_function ) 
	
      if(allocated( bundle_segment_list(i)%filter_number )) &
        DEALLOCATE( bundle_segment_list(i)%filter_number ) 
      if(allocated( bundle_segment_list(i)%Sfilter )) &
        DEALLOCATE( bundle_segment_list(i)%Sfilter ) 
      if(allocated( bundle_segment_list(i)%Zfilter )) &
        DEALLOCATE( bundle_segment_list(i)%Zfilter ) 
      if(allocated( bundle_segment_list(i)%Z_f )) &
        DEALLOCATE( bundle_segment_list(i)%Z_f ) 
      if(allocated( bundle_segment_list(i)%Zfilter_data )) &
        DEALLOCATE( bundle_segment_list(i)%Zfilter_data ) 
	    
    end do
  
    DEALLOCATE( bundle_segment_list )
  
  end if ! bundle_segment_list

  
  if (allocated( bundle_segment_geometry_list )) then
  
    do i=1,n_bundle_segment_geometries
    
      if(allocated( bundle_segment_geometry_list(i)%cable_list )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%cable_list ) 
      if(allocated( bundle_segment_geometry_list(i)%xc )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%xc ) 
      if(allocated( bundle_segment_geometry_list(i)%yc )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%yc ) 
      if(allocated( bundle_segment_geometry_list(i)%rc )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%rc ) 
      if(allocated( bundle_segment_geometry_list(i)%ri )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%ri ) 
      if(allocated( bundle_segment_geometry_list(i)%L )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%L ) 
      if(allocated( bundle_segment_geometry_list(i)%C )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%C ) 
      if(allocated( bundle_segment_geometry_list(i)%R )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%R ) 
      if(allocated( bundle_segment_geometry_list(i)%Tv )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Tv ) 
      if(allocated( bundle_segment_geometry_list(i)%Ti )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Ti ) 
      if(allocated( bundle_segment_geometry_list(i)%SC )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%SC ) 
      if(allocated( bundle_segment_geometry_list(i)%Zlink )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Zlink ) 
      if(allocated( bundle_segment_geometry_list(i)%Ylink )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Ylink ) 
      if(allocated( bundle_segment_geometry_list(i)%ZLstub )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%ZLstub ) 
      if(allocated( bundle_segment_geometry_list(i)%Yf )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Yf ) 
	
      if(allocated( bundle_segment_geometry_list(i)%filter_number )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%filter_number ) 
      if(allocated( bundle_segment_geometry_list(i)%Sfilter )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Sfilter ) 
      if(allocated( bundle_segment_geometry_list(i)%Zfilter )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Zfilter ) 
      if(allocated( bundle_segment_geometry_list(i)%Z_f )) &
        DEALLOCATE( bundle_segment_geometry_list(i)%Z_f ) 

    end do
    
    DEALLOCATE( bundle_segment_geometry_list )
  
  end if ! bundle_segment_geometry_list 
  
  
  if (allocated( cable_output_list )) then
  
    do i=1,n_cable_outputs
          if(allocated( cable_output_list(i)%conductor_list )) DEALLOCATE( cable_output_list(i)%conductor_list ) 
    end do
  
    DEALLOCATE( cable_output_list )
  
  end if ! cable_output_list

! Parallel stuff
      
  if(allocated( zmin_segment_list_send )) DEALLOCATE ( zmin_segment_list_send )
  if(allocated( zmin_segment_list_rcv  )) DEALLOCATE ( zmin_segment_list_rcv  )
  if(allocated( zmax_segment_list_send )) DEALLOCATE ( zmax_segment_list_send )
  if(allocated( zmax_segment_list_rcv  )) DEALLOCATE ( zmax_segment_list_rcv  )


  
  CALL write_line('FINISHED: deallocate_cables',0,output_to_screen_flag)

  RETURN
  
END SUBROUTINE deallocate_cables
! 
!
! NAME
!     SUBROUTINE deallocate_excitations
!
! DESCRIPTION
!     deallocate_excitations:
!
!     deallocate all allocatable structures related to geometry definition
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE deallocate_excitations

USE TLM_general
USE TLM_excitation

IMPLICIT NONE

! local variables

  integer	:: excitation_number
  integer	:: surface

! START

  
  CALL write_line('CALLED: deallocate_excitations',0,output_to_screen_flag)

  if ( allocated(excitation_functions) ) then
  
    do excitation_number=1,n_excitation_functions
  
      if ( allocated(excitation_functions(excitation_number)%value) ) then
        DEALLOCATE( excitation_functions(excitation_number)%value )
      end if
      if ( allocated(excitation_functions(excitation_number)%value_face) ) then
        DEALLOCATE( excitation_functions(excitation_number)%value_face )
      end if
  
    end do ! next excitation function

    DEALLOCATE( excitation_functions )
    
  end if   ! allocated(excitation_functions)
  
  if ( allocated( excitation_points ) ) then 
    DEALLOCATE( excitation_points ) 
  end if

  if ( allocated( excitation_surfaces ) ) then
  
    do surface=1,n_excitation_surfaces
    
      if ( allocated( excitation_surfaces(surface)%face_list ) ) then
        DEALLOCATE( excitation_surfaces(surface)%face_list )
      end if
      
      if ( allocated( excitation_surfaces(surface)%face_excitation_field_number_list ) ) then
        DEALLOCATE( excitation_surfaces(surface)%face_excitation_field_number_list )
      end if
  
    end do ! next surface
  
    DEALLOCATE( excitation_surfaces )
  
  end if
  
  if (allocated( face_excitation_field )) then
    DEALLOCATE( face_excitation_field )
  end if
  
  if (allocated( cell_excitation_field )) then
    DEALLOCATE( cell_excitation_field )
  end if
  
  if (allocated( huygens_surface%offset )) DEALLOCATE( huygens_surface%offset )
  if (allocated( huygens_surface%cx ))     DEALLOCATE( huygens_surface%cx )
  if (allocated( huygens_surface%cy ))     DEALLOCATE( huygens_surface%cy )
  if (allocated( huygens_surface%cz ))     DEALLOCATE( huygens_surface%cz )
  if (allocated( huygens_surface%face ))   DEALLOCATE( huygens_surface%face )
  if (allocated( huygens_surface%nx ))     DEALLOCATE( huygens_surface%nx )
  if (allocated( huygens_surface%ny ))     DEALLOCATE( huygens_surface%ny )
  if (allocated( huygens_surface%nz ))     DEALLOCATE( huygens_surface%nz )
  if (allocated( huygens_surface%face_excitation_field_number_list )) 	&
     DEALLOCATE( huygens_surface%face_excitation_field_number_list )

  if (allocated( excitation_mode_list )) then
  
    do surface=1,n_excitation_modes
    
      if ( allocated( excitation_mode_list(surface)%face_list ) ) then
        DEALLOCATE( excitation_mode_list(surface)%face_list )
      end if
      
      if ( allocated( excitation_mode_list(surface)%face_excitation_field_number_list ) ) then
        DEALLOCATE( excitation_mode_list(surface)%face_excitation_field_number_list )
      end if
      
      if ( allocated( excitation_mode_list(surface)%mode_field ) ) then
        DEALLOCATE( excitation_mode_list(surface)%mode_field )
      end if
    
    end do
  
    DEALLOCATE( excitation_mode_list )
  end if

  
  
  CALL write_line('FINISHED: deallocate_excitations',0,output_to_screen_flag)

  RETURN
  
END SUBROUTINE deallocate_excitations
! 
!
! NAME
!     SUBROUTINE deallocate_mode_stir
!
! DESCRIPTION
!     deallocate_mode_stir:
!
!     deallocate all allocatable structures related to geometry definition
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/01/2013 CJS
!
!
SUBROUTINE deallocate_mode_stir

USE TLM_general
USE mode_stir

IMPLICIT NONE

! local variables

  integer	:: surface

! START

  
  CALL write_line('CALLED: deallocate_mode_stir',0,output_to_screen_flag)


  if ( allocated( mode_stir_surface_list ) ) then
  
    do surface=1,n_mode_stir_surfaces
    
      if ( allocated( mode_stir_surface_list(surface)%surface_list ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%surface_list )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%surface_orientation_list ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%surface_orientation_list )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%face_list ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%face_list )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%port1 ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%port1 )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%port2 ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%port2 )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%port_voltage_list ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%port_voltage_list )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%number_of_ports_rank ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%number_of_ports_rank )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%mode_stir_voltage_pairing ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%mode_stir_voltage_pairing )
      end if
      
      if ( allocated( mode_stir_surface_list(surface)%sign ) ) then
        DEALLOCATE( mode_stir_surface_list(surface)%sign )
      end if
  
    end do ! next surface
  
    DEALLOCATE( mode_stir_surface_list )
  
  end if

  CALL write_line('FINISHED: deallocate_mode_stir',0,output_to_screen_flag)

  RETURN
  
END SUBROUTINE deallocate_mode_stir
! 
!
! NAME
!     SUBROUTINE deallocate_outputs
!
! DESCRIPTION
!     deallocate_outputs:
!
!     deallocate all allocatable structures related to geometry definition
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE deallocate_outputs

USE TLM_general
USE TLM_output

IMPLICIT NONE

! local variables

integer	:: surface
integer	:: volume

! START

  
  CALL write_line('CALLED: deallocate_outputs',0,output_to_screen_flag)
  
  if (allocated(output_points)) then
    DEALLOCATE( output_points )
  end if
  
  if (allocated(output_surfaces)) then
  
    do surface=1,n_output_surfaces
    
      if (allocated(output_surfaces(surface)%face_list)) then
        DEALLOCATE( output_surfaces(surface)%face_list )
      end if
    
      if (allocated(output_surfaces(surface)%value)) then
        DEALLOCATE( output_surfaces(surface)%value )
      end if
      
      if (allocated(output_surfaces(surface)%face_output_field_number_list)) then
        DEALLOCATE( output_surfaces(surface)%face_output_field_number_list )
      end if  
      
    end do ! next surface
    
    DEALLOCATE( output_surfaces )
    
  end if
  
  if (allocated(output_volumes)) then
  
    do volume=1,n_output_volumes
    
      if (allocated(output_volumes(volume)%cell_list)) then
        DEALLOCATE( output_volumes(volume)%cell_list )
      end if
    
      if (allocated(output_volumes(volume)%value)) then
        DEALLOCATE( output_volumes(volume)%value )
      end if
      
      if (allocated(output_volumes(volume)%cell_output_field_number_list)) then
        DEALLOCATE( output_volumes(volume)%cell_output_field_number_list )
      end if  
      
    end do ! next volume
    
    DEALLOCATE( output_volumes )
    
  end if
  
  if (allocated(output_volume_averages)) then
  
    do volume=1,n_output_volumes
    
      if (allocated(output_volume_averages(volume)%cell_list)) then
        DEALLOCATE( output_volume_averages(volume)%cell_list )
      end if
      
      if (allocated(output_volume_averages(volume)%cell_output_field_number_list)) then
        DEALLOCATE( output_volume_averages(volume)%cell_output_field_number_list )
      end if  
      
    end do ! next volume
    
    DEALLOCATE( output_volume_averages )
    
  end if
  
  if (allocated(cell_output_field)) then
    DEALLOCATE( cell_output_field )
  end if

  if (allocated(face_output_field)) then
    DEALLOCATE( face_output_field )
  end if

  if (n_far_field_surfaces.ne.0) then
  
    if (allocated( far_field_surface%face_list ))   	&
       DEALLOCATE( far_field_surface%face_list )
    if (allocated( far_field_surface%face_output_field_number_list ))   	&
       DEALLOCATE( far_field_surface%face_output_field_number_list )
    if (allocated( far_field_surface%J ))    DEALLOCATE( far_field_surface%J )
    if (allocated( far_field_surface%M ))    DEALLOCATE( far_field_surface%M )
  
  end if

  if (n_rcs_surfaces.ne.0) then
  
    if (allocated( rcs_surface%face_list ))   	&
       DEALLOCATE( rcs_surface%face_list )
    if (allocated( rcs_surface%face_output_field_number_list ))   	&
       DEALLOCATE( rcs_surface%face_output_field_number_list )
    if (allocated( rcs_surface%Etheta ))  DEALLOCATE( rcs_surface%Etheta )
    if (allocated( rcs_surface%Ephi ))    DEALLOCATE( rcs_surface%Ephi )
    if (allocated( rcs_surface%Htheta ))  DEALLOCATE( rcs_surface%Htheta )
    if (allocated( rcs_surface%Hphi ))    DEALLOCATE( rcs_surface%Hphi )
  
  end if
  
  if (n_frequency_output_surfaces.ne.0) then
  
    do surface=1,n_frequency_output_surfaces
    
      if (allocated( frequency_output_surface(surface)%face_list ))    &
         DEALLOCATE( frequency_output_surface(surface)%face_list )
      if (allocated( frequency_output_surface(surface)%face_output_field_number_list ))        &
         DEALLOCATE( frequency_output_surface(surface)%face_output_field_number_list )
      if (allocated( frequency_output_surface(surface)%value ))	&
         DEALLOCATE( frequency_output_surface(surface)%value )
  
    end do
    
    if (allocated( frequency_output_surface )) DEALLOCATE( frequency_output_surface )
  
  end if
  
  if (n_frequency_domain_power_surfaces.ne.0) then
  
    do surface=1,n_frequency_domain_power_surfaces
    
      if (allocated( frequency_domain_power_surface(surface)%face_list ))    &
         DEALLOCATE( frequency_domain_power_surface(surface)%face_list )
      if (allocated( frequency_domain_power_surface(surface)%face_output_field_number_list ))        &
         DEALLOCATE( frequency_domain_power_surface(surface)%face_output_field_number_list )
      if (allocated( frequency_domain_power_surface(surface)%E1 ))	&
         DEALLOCATE( frequency_domain_power_surface(surface)%E1 )
      if (allocated( frequency_domain_power_surface(surface)%E2 ))	&
         DEALLOCATE( frequency_domain_power_surface(surface)%E2 )
      if (allocated( frequency_domain_power_surface(surface)%H1 ))	&
         DEALLOCATE( frequency_domain_power_surface(surface)%H1 )
      if (allocated( frequency_domain_power_surface(surface)%H2 ))	&
         DEALLOCATE( frequency_domain_power_surface(surface)%H2 )
      if (allocated( frequency_domain_power_surface(surface)%Power ))	&
         DEALLOCATE( frequency_domain_power_surface(surface)%Power )
  
    end do
    
    if (allocated( frequency_domain_power_surface )) DEALLOCATE( frequency_domain_power_surface )
  
  end if
  
  if (n_SAR_volumes.ne.0) then
  
    do volume=1,n_SAR_volumes
    
      if (allocated( SAR_volume_list(volume)%cell_list ))    &
         DEALLOCATE( SAR_volume_list(volume)%cell_list )
    
      if (allocated( SAR_volume_list(volume)%Ex ))    &
         DEALLOCATE( SAR_volume_list(volume)%Ex )
    
      if (allocated( SAR_volume_list(volume)%Ey ))    &
         DEALLOCATE( SAR_volume_list(volume)%Ey )
    
      if (allocated( SAR_volume_list(volume)%Ez ))    &
         DEALLOCATE( SAR_volume_list(volume)%Ez )
  
    end do
    
    if (allocated( SAR_volume_list )) DEALLOCATE( SAR_volume_list )
  
  end if ! n_SAR_volumes.ne.0

  if (allocated( output_mode_list )) then
  
    do surface=1,n_output_modes
    
      if ( allocated( output_mode_list(surface)%face_list ) ) then
        DEALLOCATE( output_mode_list(surface)%face_list )
      end if
      
      if ( allocated( output_mode_list(surface)%face_output_field_number_list ) ) then
        DEALLOCATE( output_mode_list(surface)%face_output_field_number_list )
      end if
      
      if ( allocated( output_mode_list(surface)%mode_field ) ) then
        DEALLOCATE( output_mode_list(surface)%mode_field )
      end if
    
    end do
  
    DEALLOCATE( output_mode_list )
  end if
  
  CALL write_line('FINISHED: deallocate_outputs',0,output_to_screen_flag)

  RETURN
  
END SUBROUTINE deallocate_outputs
!
! NAME
!     deallocate_mesh
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE deallocate_mesh

USE TLM_general
USE mesh
USE TLM_excitation

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: deallocate_mesh',0,output_to_screen_flag)

  if (allocated( cell_rank )) DEALLOCATE ( cell_rank )  

  if (allocated( cell_face_rank )) DEALLOCATE ( cell_face_rank )  

  if (allocated( V )) DEALLOCATE ( V )  
  
  if (allocated( Vi_zmin )) DEALLOCATE (Vi_zmin)
  if (allocated( Vi_zmin )) DEALLOCATE (Vi_zmax)
  if (allocated( Vi_zmin )) DEALLOCATE (Vr_zmin)
  if (allocated( Vi_zmin )) DEALLOCATE (Vr_zmax)

  if (allocated( cell_centre_update_code )) DEALLOCATE ( cell_centre_update_code )

  if (allocated( face_update_code )) DEALLOCATE ( face_update_code )
  
  if (allocated( cell_update_code_to_material_data )) DEALLOCATE ( cell_update_code_to_material_data )
  if (allocated( cell_update_code_to_cable_cell_number )) DEALLOCATE ( cell_update_code_to_cable_cell_number )
  if (allocated( cell_update_code_to_excitation_number )) DEALLOCATE ( cell_update_code_to_excitation_number )
  if (allocated( cell_update_code_to_output_number )) DEALLOCATE ( cell_update_code_to_output_number )

  if (allocated( face_update_code_to_material_data )) DEALLOCATE ( face_update_code_to_material_data )
  if (allocated( face_update_code_to_cable_number )) DEALLOCATE ( face_update_code_to_cable_number )
  if (allocated( face_update_code_to_excitation_number )) DEALLOCATE ( face_update_code_to_excitation_number )
  if (allocated( face_update_code_to_output_number )) DEALLOCATE ( face_update_code_to_output_number )

  if (allocated( Vx_wrap_zmin_send )) DEALLOCATE (Vx_wrap_zmin_send)
  if (allocated( Vy_wrap_zmin_send )) DEALLOCATE (Vy_wrap_zmin_send)
  if (allocated( Vx_wrap_zmax_send )) DEALLOCATE (Vx_wrap_zmax_send)
  if (allocated( Vy_wrap_zmax_send )) DEALLOCATE (Vy_wrap_zmax_send)

  if (allocated( Vx_wrap_zmin_rcv )) DEALLOCATE (Vx_wrap_zmin_rcv)
  if (allocated( Vy_wrap_zmin_rcv )) DEALLOCATE (Vy_wrap_zmin_rcv)
  if (allocated( Vx_wrap_zmax_rcv )) DEALLOCATE (Vx_wrap_zmax_rcv)
  if (allocated( Vy_wrap_zmax_rcv )) DEALLOCATE (Vy_wrap_zmax_rcv)
  
  CALL write_line('FINISHED: deallocate_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE deallocate_mesh
