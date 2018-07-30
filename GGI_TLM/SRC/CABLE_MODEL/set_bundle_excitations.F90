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
!     23/2/2017  CJS correct the addition of source resistance and voltage for coaxial cables taking proper account
!                of the domain decomposition matrices. Before this didn't work for resitances on shield conductors.
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
  integer	:: cable_geometry_number
  integer	:: cable_type_number
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
  real*8	:: Rs1,Rs2
  integer	:: function_number_1,function_number_2
  integer       :: e1,e2

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
      
      cable_geometry_number=cable_list(cable_number)%cable_geometry_number
      cable_type_number=cable_geometry_list(cable_geometry_number)%cable_geometry_type
      
      if (rank.eq.0) then	  
        write(cable_info_file_unit,*)'Cable number',cable_number,' cable geometry number=',cable_geometry_number,  &
                                     ' cable type=',cable_type_number
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

      if (    (cable_type_number.EQ.cable_geometry_type_coaxial)           &
          .OR.(cable_type_number.EQ.cable_geometry_type_FD_coaxial) ) then

! get the elements of the matrix we are dealing with
        e1=first_conductor+1
        e2=first_conductor+2

! set source function numbers in bundle_segment data structure	
        function_number_1=cable_junction_list(junction)%excitation_function(cable_loop)%value(1)
        function_number_2=cable_junction_list(junction)%excitation_function(cable_loop)%value(2)
        
        if (rank.eq.0) then	
          write(cable_info_file_unit,*)'Coax cable source functions are:',function_number_1,function_number_1
	end if

! Check for the case where we have two different sources on inner wire and shield - the process can't deal with this at the moment
! as the sources have to be combined.
        if ( (function_number_1.NE.0).AND.(function_number_2.NE.0) ) then
        
          if (function_number_1.EQ.function_number_2) then
        
	    bundle_segment_list(segment_number)%excitation_function(e1)=0
	    bundle_segment_list(segment_number)%excitation_function(e2)=function_number_2
        
            if (rank.eq.0) then	
              write(cable_info_file_unit,*)'Setting voltage source functions'
              write(cable_info_file_unit,*)e1,bundle_segment_list(segment_number)%excitation_function(e1)
              write(cable_info_file_unit,*)e2,bundle_segment_list(segment_number)%excitation_function(e2)
	    end if
          
          else
            write(*,*)'ERROR: we cannot have two different source functions on the inner wire and the shield of a coax cable'
            write(*,*)'Junction number=',junction
            write(*,*)'Cable number=',cable_number
            write(*,*)'Cable end=',cable_end
            STOP 1
          end if
        
        else
        
	  bundle_segment_list(segment_number)%excitation_function(e1)=function_number_1-function_number_2
	  bundle_segment_list(segment_number)%excitation_function(e2)=function_number_2
        
          if (rank.eq.0) then	
            write(cable_info_file_unit,*)'Setting voltage source functions'
            write(cable_info_file_unit,*)e1,bundle_segment_list(segment_number)%excitation_function(e1)
            write(cable_info_file_unit,*)e2,bundle_segment_list(segment_number)%excitation_function(e2)
	  end if
          
        end if
        
! set resistance values in bundle_segment data structure	
        Rs1=cable_junction_list(junction)%resistance(cable_loop)%value(1)   ! value on wire
        Rs2=cable_junction_list(junction)%resistance(cable_loop)%value(2)   ! value on shield
        
! We have to take into account the domain decomposition for coaxial cables here when we add the resistances        
	bundle_segment_list(segment_number)%R(e1,e1)=bundle_segment_list(segment_number)%R(e1,e1)+Rs1+Rs2
	bundle_segment_list(segment_number)%R(e1,e2)=bundle_segment_list(segment_number)%R(e1,e2)-Rs2
	bundle_segment_list(segment_number)%R(e2,e1)=bundle_segment_list(segment_number)%R(e2,e1)-Rs2
	bundle_segment_list(segment_number)%R(e2,e2)=bundle_segment_list(segment_number)%R(e2,e2)+Rs2
        
        if (rank.eq.0) then	
          write(cable_info_file_unit,*)'Coax cable resistances are:',Rs1,Rs2
	  write(cable_info_file_unit,*)'Setting R:'
	  write(cable_info_file_unit,*)e1,e1,bundle_segment_list(segment_number)%R(e1,e1)
 	  write(cable_info_file_unit,*)e1,e2,bundle_segment_list(segment_number)%R(e1,e2)
	  write(cable_info_file_unit,*)e2,e1,bundle_segment_list(segment_number)%R(e2,e1)
	  write(cable_info_file_unit,*)e2,e2,bundle_segment_list(segment_number)%R(e2,e2)
         
	end if

      else  ! not a coaxial cable

        do conductor=1,n_conductors
      
! set source function number in bundle_segment data structure	
          excitation_function_number=cable_junction_list(junction)%excitation_function(cable_loop)%value(conductor)	
	  bundle_segment_list(segment_number)%excitation_function(first_conductor+conductor)=excitation_function_number
          
! set resistance value in bundle_segment data structure	
          Rsource=cable_junction_list(junction)%resistance(cable_loop)%value(conductor)

	  bundle_segment_list(segment_number)%R(first_conductor+conductor,first_conductor+conductor)=	&
	      bundle_segment_list(segment_number)%R(first_conductor+conductor,first_conductor+conductor)+Rsource

          if (rank.eq.0) then	
	    write(cable_info_file_unit,*)'Setting source function number ',excitation_function_number,	&
	                                 ' on conductor',first_conductor+conductor
	    write(cable_info_file_unit,*)'Setting Resistance ',Rsource,	&
	                                 ' on conductor',first_conductor+conductor
	  end if
	
        end do ! next conductor in this cable 
      
      end if   ! cable geometry type
       
    end do ! next cable loop
  
  end do ! next junction
  
  
  CALL write_line('FINISHED: set_bundle_excitations',0,output_to_screen_flag)
    
  RETURN
    
END SUBROUTINE set_bundle_excitations
