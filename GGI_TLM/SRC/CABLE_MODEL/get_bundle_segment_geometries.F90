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
!SUBROUTINE get_bundle_segment_geometries
!
! NAME
!     SUBROUTINE get_bundle_segment_geometries
!
! DESCRIPTION
!       go through all the bundle segments and work out
!       how many different geoometries there are 
!       then create and start to fill the bundle segment geometry structure
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/11/2012 CJS
!
!
SUBROUTINE get_bundle_segment_geometries()

USE TLM_general
USE Cables
USE constants
USE cell_parameters
USE File_information

IMPLICIT NONE

! local variables

  integer		:: max_cables

  integer,allocatable	:: number_of_cables_in_list(:)
  integer,allocatable	:: number_of_cables_in_list_temp(:)
  integer,allocatable	:: local_cable_list(:,:)
  integer,allocatable	:: local_cable_list_temp(:,:)
  
  integer		:: n_list1
  integer,allocatable	:: list1(:)
  integer		:: n_list2
  integer,allocatable	:: list2(:)
  
  integer		:: bundle_segment_count
  
  integer		:: n_check
  integer		:: cable
  integer		:: bundle_segment_geometry_number
  
  logical		:: cable_list_found
  
! START

  CALL write_line('CALLED: get_bundle_segment_geometries',0,output_to_screen_flag)
  
  if (n_bundle_segments.lt.1) RETURN
  
  write(cable_info_file_unit,*)
  write(cable_info_file_unit,*)'Get_bundle_segment_geometries'
  write(cable_info_file_unit,*)
  
! find the maximum number of cables in a segment - this gives one of the dimensions
! for the local_cable_list

  max_cables=0
  do bundle_segment_count=1,n_bundle_segments
  
    max_cables=max(max_cables,bundle_segment_list(bundle_segment_count)%n_cables)
    
  end do
  
  ALLOCATE( number_of_cables_in_list(1:1) )
  ALLOCATE( local_cable_list(1:1,1:max_cables) )
  ALLOCATE( list1(1:max_cables) )
  ALLOCATE( list2(1:max_cables) )
  
  local_cable_list(1:1,1:max_cables)=0
  
! set the local_cable_list to the geometry of the first bundle segment geometry

  n_bundle_segment_geometries=1
  number_of_cables_in_list(1)=bundle_segment_list(1)%n_cables
  local_cable_list( 1,1:number_of_cables_in_list(1) )=	&
        bundle_segment_list(1)%cable_list( 1:number_of_cables_in_list(1) )
  
! loop over all the bundle segments  

  do bundle_segment_count=1,n_bundle_segments

! compare the current cable list against the cable lists already found

    cable_list_found=.FALSE.

! list 1 contains the cable list which we are trying to find   
 
    n_list1=bundle_segment_list(bundle_segment_count)%n_cables
    list1(1:n_list1)=bundle_segment_list(bundle_segment_count)%cable_list(1:n_list1)
    
    do n_check=1,n_bundle_segment_geometries

! list 2 contains the cable list which we are checking against  

      n_list2=number_of_cables_in_list(n_check)
      list2(1:n_list2)=local_cable_list(n_check,1:n_list2)
      
      if (n_list1.eq.n_list2) then
! we only need to compare the lists if the number of cables is the same in each

        do cable=1,n_list1
	  if(list1(cable).ne.list2(cable)) GOTO 10
	end do

! we only get here if the cable lists are the same	
        cable_list_found=.TRUE.
	bundle_segment_geometry_number=n_check
        GOTO 20

10      CONTINUE      

      end if
 
    end do ! next list to check
    
20  if (cable_list_found) then
! set the bundle segment geometry to the geometry found
      bundle_segment_list(bundle_segment_count)%bundle_segment_geometry=bundle_segment_geometry_number
      
    else
! we need to create a new geometry in the list    

      ALLOCATE( local_cable_list_TEMP(1:n_bundle_segment_geometries,1:max_cables) )
      ALLOCATE( number_of_cables_in_list_TEMP(1:n_bundle_segment_geometries) )
      
      local_cable_list_TEMP(1:n_bundle_segment_geometries,1:max_cables)=	&
           local_cable_list(1:n_bundle_segment_geometries,1:max_cables)
	   
      number_of_cables_in_list_TEMP(1:n_bundle_segment_geometries)=	&
           number_of_cables_in_list(1:n_bundle_segment_geometries)
      
      n_bundle_segment_geometries=n_bundle_segment_geometries+1
      
      DEALLOCATE( local_cable_list )
      ALLOCATE( local_cable_list(1:n_bundle_segment_geometries,1:max_cables) )
      
      DEALLOCATE( number_of_cables_in_list )
      ALLOCATE( number_of_cables_in_list(1:n_bundle_segment_geometries) )
      
      local_cable_list(1:n_bundle_segment_geometries-1,1:max_cables)=	&
           local_cable_list_TEMP(1:n_bundle_segment_geometries-1,1:max_cables)
	   
      number_of_cables_in_list(1:n_bundle_segment_geometries-1)=	&
           number_of_cables_in_list_TEMP(1:n_bundle_segment_geometries-1)
      
      DEALLOCATE( local_cable_list_TEMP )
      DEALLOCATE( number_of_cables_in_list_TEMP )
      
      number_of_cables_in_list(n_bundle_segment_geometries)=n_list1
      local_cable_list(n_bundle_segment_geometries,1:max_cables)=0
      local_cable_list(n_bundle_segment_geometries,1:n_list1)=list1(1:n_list1)
      
! set the bundle segment geometry to the new geometry 
      bundle_segment_geometry_number=n_bundle_segment_geometries
      bundle_segment_list(bundle_segment_count)%bundle_segment_geometry=bundle_segment_geometry_number
      
    end if    
    
  end do ! next bundle segment requiring a bundle geometry 

! create the  bundle_segment_geometry structure

  ALLOCATE( bundle_segment_geometry_list(1:n_bundle_segment_geometries) )
  
  write(cable_info_file_unit,*)'Number of bundle segment geometries=',n_bundle_segment_geometries
  
  do bundle_segment_geometry_number=1,n_bundle_segment_geometries
  
     n_list1=number_of_cables_in_list(bundle_segment_geometry_number)
     
     bundle_segment_geometry_list(bundle_segment_geometry_number)%n_cables=n_list1
     
     ALLOCATE( bundle_segment_geometry_list(bundle_segment_geometry_number)%cable_list(1:n_list1) )
     
     bundle_segment_geometry_list(bundle_segment_geometry_number)%cable_list(1:n_list1)=	&
           local_cable_list(bundle_segment_geometry_number,1:n_list1)
	   
     write(cable_info_file_unit,*)'Bundle segment geometry number=',bundle_segment_geometry_number
     write(cable_info_file_unit,*)'Number of cables=',	&
           bundle_segment_geometry_list(bundle_segment_geometry_number)%n_cables
     write(cable_info_file_unit,*)'Cable list:',	&
           bundle_segment_geometry_list(bundle_segment_geometry_number)%cable_list(1:n_list1)
	   
  end do ! next bundle_segment_geometry
  
  DEALLOCATE( number_of_cables_in_list )
  DEALLOCATE( local_cable_list )
  DEALLOCATE( list1 )
  DEALLOCATE( list2 )

  CALL write_line('FINISHED: get_bundle_segment_geometries',0,output_to_screen_flag)
    
  RETURN
   
9000 CALL write_line('Error in get_bundle_segment_geometries:',0,.TRUE.)
     CALL write_line('',0,.TRUE.)
     STOP
  
  
END SUBROUTINE get_bundle_segment_geometries
