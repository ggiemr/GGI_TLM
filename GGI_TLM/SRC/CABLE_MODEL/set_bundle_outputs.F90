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
!SUBROUTINE set_bundle_outputs
!
! NAME
!     SUBROUTINE set_bundle_outputs
!
! DESCRIPTION
!
!     
! COMMENTS
!     
!
!  
!
! HISTORY
!
!     started 21/09/2012 CJS
!
!
SUBROUTINE set_bundle_outputs()

USE TLM_general
USE geometry_types
USE Cables
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer		:: output_number
  integer		:: closest_point_number
  integer		:: cable_number
  type(cell_point)	:: output_point
  
  integer		:: segment_number
  
  integer		:: cable_loop
  integer		:: cable
  integer		:: n_cables_found
  integer		:: cable_geometry
  integer		:: conductor
  integer		:: n_conductors
  integer		:: bundle_segment_conductor
  integer		:: output_conductor_count

! START

  CALL write_line('CALLED: set_bundle_outputs',0,output_to_screen_flag)
  
  write(cable_info_file_unit,*)'Set bundle outputs'

! loop over cable outputs
  do output_number=1,n_cable_outputs
    
    write(cable_info_file_unit,*)'Output number:',output_number
  
    closest_point_number=cable_output_list(output_number)%closest_point_number
    cable_number=cable_output_list(output_number)%cable_number
    
! find the closest bundle segment containing the required cable_number
    CALL closest_cable_segment(closest_point_number,cable_number,segment_number,output_point)
    
    write(cable_info_file_unit,*)'Found closest bundle segment:',segment_number
    write(cable_info_file_unit,*)'Cell_point                  :',output_point

    cable_output_list(output_number)%bundle_segment_number=segment_number
    cable_output_list(output_number)%output_point=output_point
    
    n_cables_found=0
    do cable=1,bundle_segment_list(segment_number)%n_cables
      if (bundle_segment_list(segment_number)%cable_list(cable).eq.cable_number) n_cables_found=n_cables_found+1
    end do
    
    write(cable_info_file_unit,*)'Number of cables found:',n_cables_found
    
    cable_geometry=cable_list(cable_number)%cable_geometry_number
    
    cable_output_list(output_number)%n_conductors=n_cables_found*cable_geometry_list(cable_geometry)%n_conductors  

! loop through all the cables in the output bundle_segment and build the cable_output conductor_list
    ALLOCATE( cable_output_list(output_number)%conductor_list(1:cable_output_list(output_number)%n_conductors) )
     
    n_cables_found=0
    output_conductor_count=0
    bundle_segment_conductor=0
    
    do cable_loop=1,bundle_segment_list(segment_number)%n_cables
    
      cable=bundle_segment_list(segment_number)%cable_list(cable_loop)
      cable_geometry=cable_list(cable)%cable_geometry_number
      n_conductors=cable_geometry_list(cable_geometry)%n_conductors
      
      if (cable.eq.cable_number) then
      
        n_cables_found=n_cables_found+1
	do conductor=1,n_conductors
	  output_conductor_count=output_conductor_count+1
          bundle_segment_conductor=bundle_segment_conductor+1
	  cable_output_list(output_number)%conductor_list(output_conductor_count)=bundle_segment_conductor
	end do
	
      else
      
        bundle_segment_conductor=bundle_segment_conductor+n_conductors
        
      end if ! is this the required cable
      
    end do
    
! counting checks
    if ( bundle_segment_conductor.NE.bundle_segment_list(segment_number)%n_conductors ) GOTO 9000
    if ( output_conductor_count.NE.cable_output_list(output_number)%n_conductors )      GOTO 9010
    
  end do ! next cable output

  CALL write_line('FINISHED: set_bundle_outputs',0,output_to_screen_flag)
    
  RETURN
  
9000 CALL write_line('ERROR in set bundle outputs',0,.TRUE.)
     CALL write_line('Mismatch in bundle segment conductor count',0,.TRUE.)
     CALL write_line_integer('bundle_segment_list(segment_number)%n_conductors',	&
                              bundle_segment_list(segment_number)%n_conductors,0,.TRUE.)
     CALL write_line_integer('bundle_segment_conductor count',bundle_segment_conductor,0,.TRUE.)
     STOP
     
9010 CALL write_line('ERROR in set bundle outputs',0,.TRUE.)
     CALL write_line('Mismatch in output conductor count',0,.TRUE.)
     CALL write_line_integer('cable_output_list(output_number)%n_conductors',	&
                              cable_output_list(output_number)%n_conductors,0,.TRUE.)
     CALL write_line_integer('output_conductor_count',output_conductor_count,0,.TRUE.)
     STOP
  
END SUBROUTINE set_bundle_outputs
