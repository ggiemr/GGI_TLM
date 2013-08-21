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
! SUBROUTINE check_cable_input_data
!
! NAME
!     check_cable_input_data
!
! DESCRIPTION
!     reset all packet data
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 12/10/2012 CJS
!
!
SUBROUTINE check_cable_input_data

USE TLM_general
USE geometry
USE Cables

IMPLICIT NONE

! local variables

  integer	:: i,j,k
  integer	:: line
  integer	:: cable
  integer	:: cable_end
  integer	:: cable_geometry
  
  integer 	:: n_conductors
  integer 	:: n_conductors2

! START  

  CALL write_line('CALLED: check_cable_input_data',0,output_to_screen_flag)

! check that the cable geometry exists and also whether all the lines on cable routes exist.  
  do i=1,n_cables
  
    cable_geometry=cable_list(i)%cable_geometry_number
    if ( (cable_geometry.gt.n_cable_geometries).OR.(cable_geometry.lt.1) )then
    
      CALL write_line('ERROR in cable list: cable geometry does not exist',0,.TRUE.)
      CALL write_line_integer('cable number:',i,0,.TRUE.)
      CALL write_line_integer('cable geometry number:',cable_geometry,0,.TRUE.)
      CALL write_line_integer('Number of cable geometries=',n_cable_geometries,0,.TRUE.)
      STOP

    end if
    
    do j=1,cable_list(i)%n_lines
    
      line=cable_list(i)%line_list(j)
      if ( (line.gt.n_lines).OR.(line.lt.1) )then
    
        CALL write_line('ERROR in cable list: line does not exist',0,.TRUE.)
        CALL write_line_integer('cable number:',i,0,.TRUE.)
        CALL write_line_integer('line number:',line,0,.TRUE.)
        CALL write_line_integer('Number of lines=',n_lines,0,.TRUE.)
        STOP

      end if
      
    end do ! next line on cable route
    
  end do ! next cable

! check whether the cable required for output exists  
  do i=1,n_cable_outputs
    cable=cable_output_list(i)%cable_number
    if ( (cable.gt.n_cables).OR.(cable.lt.1) )then
    
      CALL write_line('ERROR in cable_output_list: cable does not exist',0,.TRUE.)
      CALL write_line_integer('cable output number:',i,0,.TRUE.)
      CALL write_line_integer('cable number:',cable,0,.TRUE.)
      CALL write_line_integer('Number of cables=',n_cables,0,.TRUE.)
      STOP

    end if
  end do

! check that cables at a junction exist, the end number is 1 or 2 
! also check the dimensions of the Pmatrices specified
  do i=1,n_cable_junctions
      
    do j=1,cable_junction_list(i)%number_of_cables
    
      cable=cable_junction_list(i)%cable_list(j)
      if ( (cable.gt.n_cables).OR.(cable.lt.1) )then
    
        CALL write_line('ERROR in cable_junction_list: cable does not exist',0,.TRUE.)
        CALL write_line_integer('cable junction number:',i,0,.TRUE.)
        CALL write_line_integer('cable number:',cable,0,.TRUE.)
        CALL write_line_integer('Number of cables=',n_cables,0,.TRUE.)
        STOP

      end if
    
      cable_end=cable_junction_list(i)%cable_end_list(j)
      if ( (cable_end.ne.1).AND.(cable_end.ne.2) )then
    
        CALL write_line('ERROR in cable_junction_list: cable end should be 1 or 2',0,.TRUE.)
        CALL write_line_integer('cable junction number:',i,0,.TRUE.)
        CALL write_line_integer('cable number:',cable,0,.TRUE.)
        CALL write_line_integer('cable end:',cable_end,0,.TRUE.)
        STOP

      end if

! n_conductors as set in the junction information      
      n_conductors=cable_junction_list(i)%n_external_conductors(j)
      cable_geometry=cable_list(cable)%cable_geometry_number
      n_conductors2=cable_geometry_list(cable_geometry)%n_conductors
      
      if (n_conductors.ne.n_conductors2) then
    
        CALL write_line('ERROR in cable_junction_list: Discrepancy in number of conductors',0,.TRUE.)
        CALL write_line_integer('cable junction number:',i,0,.TRUE.)
        CALL write_line_integer('cable number:',cable,0,.TRUE.)
        STOP

      end if
      
      
    end do ! next cable in cable list
    
  end do ! next cable junction

  CALL write_line('FINISHED: check_cable_input_data',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE check_cable_input_data
