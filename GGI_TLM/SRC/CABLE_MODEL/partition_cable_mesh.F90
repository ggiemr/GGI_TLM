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
!
! Name partition_cable_mesh
!     
!
! Description
!     
!
! Comments:
!      
!
! History
!
!     started 12/05/09 CJS
!     Adapted for GGI_TLM from Fieldsolve code 27/11/2012 CJS
!     Revised for more efficient two stage parallel process 13/12/2012 CJS
!

  SUBROUTINE partition_cable_mesh()

USE TLM_general
USE mesh
USE cell_parameters
USE cables
USE file_information

IMPLICIT NONE

! variables passed to subroutine


! local_variables

  integer :: cz
  integer :: zmin_count_send,zmax_count_send
  integer :: zmin_count_rcv,zmax_count_rcv
  integer :: output
  integer :: i
  
  integer :: face_segment
  integer :: cell_junction
  
! function_types

! START
	   
  CALL write_line('FINISHED: partition_cable_mesh',0,timestepping_output_to_screen_flag)
 
  n_zmin_segments_send=0
  n_zmax_segments_send=0
  n_zmin_segments_rcv=0
  n_zmax_segments_rcv=0
 
  do cell_junction=1,n_cell_centre_junctions
  
! get the TLM cell z coordinate and count interface cells
    cz=cell_centre_junction_list(cell_junction)%cell_point%cell%k
	  
    if ( (cz.eq.nz1-1).AND.(rank.NE.0) ) then       
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are received from the rank-1 process
        n_zmin_segments_rcv=n_zmin_segments_rcv+1
      end if          
    end if
        
    if ( (cz.eq.nz1-1).AND.(rank.NE.0) ) then       
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are sent to the rank-1 process
        n_zmin_segments_send=n_zmin_segments_send+1
      end if          
    end if
        
    if ( (cz.eq.nz2).AND.(rank.NE.np-1) ) then        
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are sent to the rank+1 process
        n_zmax_segments_send=n_zmax_segments_send+1
      end if          
    end if
        
    if ( (cz.eq.nz2).AND.(rank.NE.np-1) ) then        
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are received from the rank+1 process
        n_zmax_segments_rcv=n_zmax_segments_rcv+1
      end if          
    end if
	
  end do
      
  ALLOCATE ( zmin_segment_list_send(1:n_zmin_segments_send) )
  ALLOCATE ( zmin_segment_list_rcv (1:n_zmin_segments_rcv) )
  ALLOCATE ( zmax_segment_list_send(1:n_zmax_segments_send) )
  ALLOCATE ( zmax_segment_list_rcv (1:n_zmax_segments_rcv) )
  
  zmin_count_send=0
  zmin_count_rcv=0
  zmax_count_send=0
  zmax_count_rcv=0
  n_zmin_reals_send=0
  n_zmin_reals_rcv=0
  n_zmax_reals_send=0
  n_zmax_reals_rcv=0

! generate list of segments which lie on the interface between processes
! for sending and receiveing
 
  do cell_junction=1,n_cell_centre_junctions
  
! get the TLM cell z coordinate       
    cz=cell_centre_junction_list(cell_junction)%cell_point%cell%k
	  
    if ( (cz.eq.nz1-1).AND.(rank.NE.0) ) then       
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are received from the rank-1 process
    	zmin_count_rcv=zmin_count_rcv+1
    	zmin_segment_list_rcv(zmin_count_rcv)=face_segment
    	n_zmin_reals_rcv=n_zmin_reals_rcv+bundle_segment_list(face_segment)%n_conductors
      end if          
    end if
        
    if ( (cz.eq.nz1-1).AND.(rank.NE.0) ) then       
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are sent to the rank-1 process
    	zmin_count_send=zmin_count_send+1
    	zmin_segment_list_send(zmin_count_send)=face_segment
    	n_zmin_reals_send=n_zmin_reals_send+bundle_segment_list(face_segment)%n_conductors
      end if          
    end if
        
    if ( (cz.eq.nz2).AND.(rank.NE.np-1) ) then        
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are sent to the rank+1 process
    	zmax_count_send=zmax_count_send+1
    	zmax_segment_list_send(zmax_count_send)=face_segment
    	n_zmax_reals_send=n_zmax_reals_send+bundle_segment_list(face_segment)%n_conductors
      end if          
    end if
        
    if ( (cz.eq.nz2).AND.(rank.NE.np-1) ) then        
      face_segment=cell_centre_junction_list(cell_junction)%segment_list(face_zmax)
      if (face_segment.ne.0) then
! these segments are received from the rank+1 process
    	zmax_count_rcv=zmax_count_rcv+1
    	zmax_segment_list_rcv(zmax_count_rcv)=face_segment
    	n_zmax_reals_rcv=n_zmax_reals_rcv+bundle_segment_list(face_segment)%n_conductors
      end if          
    end if

  end do ! next cell_junction

  allocate(wire_Vr_zmin(1:n_zmin_reals_send))
  allocate(wire_Vi_zmin(1:n_zmin_reals_rcv))
  allocate(wire_Vr_zmax(1:n_zmax_reals_send))
  allocate(wire_Vi_zmax(1:n_zmax_reals_rcv))      

  do output=1,n_cable_outputs
  
    cz=cable_output_list(output)%output_point%cell%k
    cable_output_list(output)%rank=cell_rank(cz)
    	       
  end do  ! next cable
  	   
  CALL write_line('FINISHED: partition_cable_mesh',0,timestepping_output_to_screen_flag)
	           
  RETURN
      
END SUBROUTINE partition_cable_mesh



