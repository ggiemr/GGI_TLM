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
!SUBROUTINE create_bundle_segments
!
! NAME
!     SUBROUTINE create_bundle_segments
!
! DESCRIPTION
!       build up bundle segments from individual cable meshes
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
SUBROUTINE create_bundle_segments()

USE TLM_general
USE Cables
USE constants
USE cell_parameters
USE File_information

IMPLICIT NONE

! local variables

  integer,allocatable	:: cable_count(:,:,:,:)
  integer,allocatable	:: bundle_segment_number(:,:,:,:)
  
  integer		:: cable
  integer		:: segment
  type(cell_point)	:: cable_cell_segment
  
  integer		:: ix,iy,iz,face
  integer		:: bundle_segment_count
  integer		:: number_of_cables
  
! START

  CALL write_line('CALLED: create_bundle_segments',0,output_to_screen_flag)
  
! allocate a full mesh to count the total number of bundle segments
  ALLOCATE( cable_count(1:nx,1:ny,1:nz,1:6) )
  ALLOCATE( bundle_segment_number(1:nx,1:ny,1:nz,1:6) )

  cable_count(1:nx,1:ny,1:nz,1:6)=0
  bundle_segment_number(1:nx,1:ny,1:nz,1:6)=0

! loop over cables adding each cell_segment into the cable_bundle_segment_mesh count
  n_bundle_segments=0
  
!  write(cable_info_file_unit,*)'Setting cable mesh'
  
  do cable=1,n_cables
  
    do segment=1,cable_list(cable)%number_of_cable_segments
    
!      write(cable_info_file_unit,*)'Setting cable,',cable,' segment',segment,' in mesh'

! get the cell_point which is on the cell face  
      if (cable_list(cable)%cable_segment_list(segment)%segment_point(1)%point.NE.centre) then
        cable_cell_segment=cable_list(cable)%cable_segment_list(segment)%segment_point(1)
      else if (cable_list(cable)%cable_segment_list(segment)%segment_point(2)%point.NE.centre) then
        cable_cell_segment=cable_list(cable)%cable_segment_list(segment)%segment_point(2)
      else
        GOTO 9000 ! error, no cell face point found
      end if
      
      ix=cable_cell_segment%cell%i
      iy=cable_cell_segment%cell%j
      iz=cable_cell_segment%cell%k
      face=cable_cell_segment%point
      
      if (cable_count(ix,iy,iz,face).eq.0) then
! new bundle segment found on this face
        n_bundle_segments=n_bundle_segments+1
      end if
            
! count the new cable in this segment
      cable_count(ix,iy,iz,face)=cable_count(ix,iy,iz,face)+1
      
    end do ! next cable segment
  
  end do ! next cable
  
  write(cable_info_file_unit,*)'n_bundle_segments',n_bundle_segments
  
! Allocate the bundle_segment_list
  ALLOCATE( bundle_segment_list(1:n_bundle_segments) )
  
! loop over the mesh, allocating the cable_lists for each bundle segment

  write(cable_info_file_unit,*)'Allocating cable lists'
  bundle_segment_count=0
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx
        do face=1,6
	
          if (cable_count(ix,iy,iz,face).ne.0) then

! Allocate cable_list for this bundle segment	  
	    bundle_segment_count=bundle_segment_count+1
	    
!	    write(cable_info_file_unit,*)'Allocating bundle_segment',bundle_segment_count
	    
	    number_of_cables=cable_count(ix,iy,iz,face)
	    bundle_segment_list(bundle_segment_count)%n_cables=number_of_cables
	    ALLOCATE( bundle_segment_list(bundle_segment_count)%cable_list(1:number_of_cables) )
            bundle_segment_list(bundle_segment_count)%cable_list(1:number_of_cables)=0
	    
! set cable_bundle_segment_mesh to the bundle_segment number
            bundle_segment_number(ix,iy,iz,face)=bundle_segment_count
	    
	  end if
	  
        end do 
      end do 
    end do 
  end do 
  
! Loop over the cable list again setting the cable numbers in the bundle segment list

  cable_count(1:nx,1:ny,1:nz,1:6)=0
  
  write(cable_info_file_unit,*)'Setting cable lists'
  
  do cable=1,n_cables
  
    do segment=1,cable_list(cable)%number_of_cable_segments

! get the cell_point which is on the cell face  
      if (cable_list(cable)%cable_segment_list(segment)%segment_point(1)%point.NE.centre) then
        cable_cell_segment=cable_list(cable)%cable_segment_list(segment)%segment_point(1)
      else if (cable_list(cable)%cable_segment_list(segment)%segment_point(2)%point.NE.centre) then
        cable_cell_segment=cable_list(cable)%cable_segment_list(segment)%segment_point(2)
      else
        GOTO 9000 ! error, no cell face point found
      end if
      
      ix=cable_cell_segment%cell%i
      iy=cable_cell_segment%cell%j
      iz=cable_cell_segment%cell%k
      face=cable_cell_segment%point
	   
! count the new cable in this segment
      cable_count(ix,iy,iz,face)=cable_count(ix,iy,iz,face)+1
         
!      write(cable_info_file_unit,*)'Setting cable in bundle_segment',bundle_segment_number(ix,iy,iz,face)
!      write(cable_info_file_unit,*)'Cable=',cable

! add the new cable to the bundle segment
      bundle_segment_list(bundle_segment_number(ix,iy,iz,face))%cable_list(cable_count(ix,iy,iz,face))=cable

      bundle_segment_list(bundle_segment_number(ix,iy,iz,face))%cable_segment=	&
                          cable_list(cable)%cable_segment_list(segment)
			  
! add the bundle segment to the cable based bundle_segment list
      cable_list(cable)%bundle_segment_list(segment)=bundle_segment_number(ix,iy,iz,face)
      
    end do ! next cable segment
  
  end do ! next cable
  
!  write(cable_info_file_unit,*)''
!  write(cable_info_file_unit,*)'Bundle segment list'
!  
!  do bundle_segment_count=1,n_bundle_segments
!
!    write(cable_info_file_unit,*)bundle_segment_count,bundle_segment_list(bundle_segment_count)%cable_list(:)
!
!  end do
  
  
  DEALLOCATE( cable_count )
  DEALLOCATE( bundle_segment_number )

  CALL write_line('FINISHED: create_bundle_segments',0,output_to_screen_flag)
    
  RETURN
   
9000 CALL write_line('Error in create_bundle_segments:',0,.TRUE.)
     CALL write_line('No cell face point found',0,.TRUE.)
     STOP
  
  
END SUBROUTINE create_bundle_segments
