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
!SUBROUTINE set_bundle_excitations
!
! NAME
!     SUBROUTINE set_bundle_excitations
!
! DESCRIPTION
!
!     
! COMMENTS
!     
!
!  Set excitations into the bundle_segment data structure
!  Excitation data is defined in the cable_junction packet
!
! HISTORY
!
!     started 21/09/2012 CJS
!
!
SUBROUTINE set_bundle_excitations()

USE TLM_general
USE Cables
USE constants
USE file_information

IMPLICIT NONE

! local variables

  integer	:: junction
  integer	:: cable_loop
  integer	:: cable_number
  integer	:: segment_cable_loop
  integer	:: cable_end
  integer	:: n_cables_found
  integer	:: number_of_cable_segments
  integer	:: segment_number
  integer	:: segment_cable
  integer	:: cable_geometry_segment_cable
  integer	:: n_conductors_segment_cable
  integer	:: first_conductor
  integer	:: n_conductors,conductor
  integer	:: excitation_function_number
  real*8	:: Rsource

! START

  CALL write_line('CALLED: set_bundle_excitations',0,output_to_screen_flag)

  if (rank.eq.0) then	    
    write(cable_info_file_unit,*)'Bundle excitations'    
  end if
  
  do junction=1,n_cable_junctions
  
    if (rank.eq.0) then	  
      write(cable_info_file_unit,*)'Junction number',junction   
    end if
  
    do cable_loop=1,cable_junction_list(junction)%number_of_cables
      
      cable_number=cable_junction_list(junction)%cable_list(cable_loop)
      cable_end=cable_junction_list(junction)%cable_end_list(cable_loop)
      n_conductors=cable_junction_list(junction)%n_external_conductors(cable_loop)
  
      if (rank.eq.0) then	  
        write(cable_info_file_unit,*)'Cable number',cable_number
      end if
      
! get the bundle segment for this cable end
      if (cable_end.eq.1) then
        segment_number=cable_list(cable_number)%bundle_segment_list(1)
      else
        number_of_cable_segments=cable_list(cable_number)%number_of_cable_segments
        segment_number=cable_list(cable_number)%bundle_segment_list(number_of_cable_segments)
      end if
      
      n_cables_found=0
      do segment_cable_loop=1,bundle_segment_list(segment_number)%n_cables
        if (bundle_segment_list(segment_number)%cable_list(segment_cable_loop).eq.cable_number) n_cables_found=n_cables_found+1
      end do
      
      if (n_cables_found.NE.1) then
      
! write a warning 
        CALL write_line('ERROR in set_bundle_excitations',warning_file_unit,.TRUE.)
        CALL write_line('Cable appears in bundle_segment more than once',warning_file_unit,.TRUE.)
        CALL write_line_integer('n_cables_found=',n_cables_found,warning_file_unit,.TRUE.)

        CALL write_line('ERROR in set_bundle_excitations',cable_info_file_unit,.TRUE.)
        CALL write_line('Cable appears in bundle_segment more than once',cable_info_file_unit,.TRUE.)
        CALL write_line_integer('n_cables_found=',n_cables_found,cable_info_file_unit,.TRUE.)

        CALL write_line('ERROR in set_bundle_excitations',0,.TRUE.)
        CALL write_line('Cable appears in bundle_segment more than once',0,.TRUE.)
        CALL write_line_integer('n_cables_found=',n_cables_found,0,.TRUE.)

      end if

! find the first conductor number in the bundle_segment which relates to this cable
      first_conductor=0
      do segment_cable_loop=1,bundle_segment_list(segment_number)%n_cables
    
        segment_cable=bundle_segment_list(segment_number)%cable_list(segment_cable_loop)
        cable_geometry_segment_cable=cable_list(segment_cable)%cable_geometry_number
        n_conductors_segment_cable=cable_geometry_list(cable_geometry_segment_cable)%n_conductors
      
        if (segment_cable.eq.cable_number) GOTO 1000
      	
        first_conductor=first_conductor+n_conductors_segment_cable
      
      end do

1000 CONTINUE ! jump out of loop to here if cable is found      
      do conductor=1,n_conductors
      
        excitation_function_number=cable_junction_list(junction)%excitation_function(cable_loop)%value(conductor)
        Rsource=cable_junction_list(junction)%resistance(cable_loop)%value(conductor)
	
! set source function number in bundle_segment data structure	
	bundle_segment_list(segment_number)%excitation_function(first_conductor+conductor)=excitation_function_number

! set resistance value in bundle_segment data structure	
	bundle_segment_list(segment_number)%R(first_conductor+conductor,first_conductor+conductor)=	&
	    bundle_segment_list(segment_number)%R(first_conductor+conductor,first_conductor+conductor)+Rsource

        if (rank.eq.0) then	
	  write(cable_info_file_unit,*)'Setting source function number ',excitation_function_number,	&
	                               ' on conductor',first_conductor+conductor
	  write(cable_info_file_unit,*)'Setting Resistance ',Rsource,	&
	                               ' on conductor',first_conductor+conductor
	end if
	
      end do ! next conductor in this cable 
      
    end do ! next cable loop
  
  end do ! next junction
  
  
  CALL write_line('FINISHED: set_bundle_excitations',0,output_to_screen_flag)
    
  RETURN
    
END SUBROUTINE set_bundle_excitations
