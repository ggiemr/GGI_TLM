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
! SUBROUTINE read_mode_stir_surface_list
!
! NAME
!     read_mode_stir_surface_list
!
! DESCRIPTION
!     read mode stir surface list packet
!
! Example packet:
!
!Mode_stir_surface_list
!2   Number of mode_stir_surfaces (integer)
!1	MODE_STIR_SURFACE NUMBER
!0.999  	! boundary condition voltage reflection coefficient 
!0.0	! amplitude of random impulse excitation on the boundary
!1	! number of geometric surfaces in this mode_stir_surface
!1	! surface list
!1	! surface orientation list
!2	MODE_STIR_SURFACE NUMBER
!0.999  	! boundary condition voltage reflection coefficient 
!1.0	! amplitude of random impulse excitation on the boundary
!1	! number of geometric surfaces in this mode_stir_surface
!0	! surface list - surface 0 is the outer boundary
!1	! surface orientation list
!
! COMMENTS
!     
!
! HISTORY
!
!     started 10/01/2013 CJS
!
!
SUBROUTINE read_mode_stir_surface_list

USE TLM_general
USE mode_stir
USE file_information
USE geometry
USE constants

IMPLICIT NONE

! local variables

integer	:: surface_number
integer :: read_number

integer	:: n_geometric_surfaces
integer	:: i

! START  

  CALL write_line('CALLED: read_mode_stir_surface_list',0,output_to_screen_flag)

  read(input_file_unit,*,err=9005)n_mode_stir_surfaces
  
  CALL write_line_integer('number of mode stir surfaces',n_mode_stir_surfaces,0,output_to_screen_flag)
  
  if ( allocated( mode_stir_surface_list ) ) GOTO 9000
  
  ALLOCATE ( mode_stir_surface_list(1:n_mode_stir_surfaces) )

  do surface_number=1,n_mode_stir_surfaces
  
    CALL write_line_integer('Reading surface number',surface_number,0,output_to_screen_flag)
    
    read(input_file_unit,*,err=9005)read_number
    if (read_number.ne.surface_number) goto 9010

    read(input_file_unit,*,err=9005)mode_stir_surface_list(surface_number)%R_ms
    read(input_file_unit,*,err=9005)mode_stir_surface_list(surface_number)%impulse_amplitude
    
    read(input_file_unit,*,err=9005)mode_stir_surface_list(surface_number)%number_of_surfaces
    n_geometric_surfaces=mode_stir_surface_list(surface_number)%number_of_surfaces
    
    ALLOCATE(mode_stir_surface_list(surface_number)%surface_list(1:n_geometric_surfaces))
    ALLOCATE(mode_stir_surface_list(surface_number)%surface_orientation_list(1:n_geometric_surfaces))
    
    read(input_file_unit,*,err=9020)( mode_stir_surface_list(surface_number)%surface_list(i)	&
                                     ,i=1,n_geometric_surfaces )
    read(input_file_unit,*,err=9030)( mode_stir_surface_list(surface_number)%surface_orientation_list(i)	&
                                     ,i=1,n_geometric_surfaces )

! check orientation flag is +1 or -1				     
    do i=1,n_geometric_surfaces
      if (abs(mode_stir_surface_list(surface_number)%surface_orientation_list(i)).ne.1) goto 9040
    end do
    
  end do ! next mode_stir_surface

  CALL write_line('FINISHED: read_mode_stir_surface_list',0,output_to_screen_flag)
  
  RETURN
  
9000 CALL write_line('Error allocating mode_stir_surfaces:',0,.TRUE.)
     CALL write_line('mode_stir_surfaces already allocated',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP
  
9005 CALL write_line('Error reading mode_stir_surface_list packet from input file:',0,.TRUE.)
     CALL write_error_line(input_file_unit)
     STOP

9010 CALL write_line('Error reading mode_stir_surface_list packet data',0,.TRUE.)
     CALL write_line('Surfaces should be numbered in order at the moment...',0,.TRUE.)
     STOP
  
9020 CALL write_line('Error reading mode_stir_surface_list_packet_data',0,.TRUE.)
     CALL write_line('Problem in reading surface_list',0,.TRUE.)
     STOP
    
9030 CALL write_line('Error reading mode_stir_surface_list_packet_data',0,.TRUE.)
     CALL write_line('Problem in reading surface_orientation_list',0,.TRUE.)
     STOP
    
9040 CALL write_line('Error reading mode_stir_surface_list_packet_data',0,.TRUE.)
     CALL write_line('Surface_orientation should be +1 or -1',0,.TRUE.)
     STOP
 
END SUBROUTINE read_mode_stir_surface_list
