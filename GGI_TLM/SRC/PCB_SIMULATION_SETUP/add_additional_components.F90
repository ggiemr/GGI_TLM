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
! NAME
!     SUBROUTINE specify_additional_components
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
!     started 14/3/19 CJS
!     
!

SUBROUTINE specify_additional_components

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i

character(LEN=256) :: component_type

real :: hxmin,hxmax,hymin,hymax,hzmin,hzmax
integer :: nslots
real*8  :: ws,wr,ds

integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
real*8  :: lx,ly,lz
integer :: ix,iy,iz,slot

real*8 :: slot_xmin,slot_xmax,slot_ymin,slot_ymax,slot_zmin,slot_zmax

integer :: dslot

! START

  write(*,*)'CALLED specify_additional_components'

  read(10,*)n_additional_components
  
  do i=1,n_additional_components
  
    read(10,'(A)')component_type
    
    if( index(component_type,'heatsink').NE.0 ) then
    
      read(10,*)hxmin,hymin,hzmin,hxmax,hymax,hzmax
    
      read(10,*)nslots
      
      read(10,*)ws
      
      read(10,*)ds
    
      lx=hxmax-hxmin
      ly=hymax-hymin
      lz=hzmax-hzmin
    
! calculate the extent of the heatsinkk in cells
    
      CALL get_TLM_cell_from_coordinate(hxmin+dl/2d0,hymin+dl/2d0,hzmin+dl/2d0,ixmin,iymin,izmin)
      CALL get_TLM_cell_from_coordinate(hxmax-dl/2d0,hymax-dl/2d0,hzmax-dl/2d0,ixmax,iymax,izmax)
      
      write(*,*)'Heatsink minimum coordinates:'
      write(*,*)hxmin,hymin,hzmin
      write(*,*)'Heatsink maximum coordinates:'
      write(*,*)hxmax,hymax,hzmax
      
      write(*,*)'Heatsink extent:'
      write(*,*)lx,ly,lz
      
      write(*,*)'Heatsink minimum cell:'
      write(*,*)ixmin,iymin,izmin
      write(*,*)'Heatsink maximum cell:'
      write(*,*)ixmax,iymax,izmax

! set the whole heatsink volume to PEC
    
      do ix=ixmin,ixmax
        do iy=iymin,iymax
          do iz=izmin,izmax
            material_mesh(centre,ix,iy,iz)=1
          end do
        end do
      end do

! set the slots in the heatsink to free space

      wr=(lx-nslots*wr)/dble(nslots+1)
      if (wr.LT.0d0) then
        write(*,*)'ERROR: width of ribs in heatsink is less than zero'
        write(*,*)'wr=',wr
        STOP 1
      end if
      
      write(*,*)'Heatsink slot width:',ws
      write(*,*)'Heatsink rib width :',wr
      write(*,*)'Heatsink depth     :',ds

! assume slots are in the zmin face and are in the y direction

      slot_ymin=hymin
      slot_ymax=hymax
      slot_zmin=hzmin
      slot_zmax=hzmin+ds
      
      do slot=1,nslots
      
        slot_xmin=hxmin+wr+(slot-1)*(wr+ws)
        slot_xmax=slot_xmin+ws
        
        write(*,*)'Slot number:',slot
        write(*,*)'xmin:',slot_xmin,' xmax:',slot_xmax
        
        CALL get_TLM_cell_from_coordinate(slot_xmin+dl/2d0,slot_ymin+dl/2d0,slot_zmin+dl/2d0,ixmin,iymin,izmin)
        CALL get_TLM_cell_from_coordinate(slot_xmax-dl/2d0,slot_ymax-dl/2d0,slot_zmax-dl/2d0,ixmax,iymax,izmax)

        write(*,*)'From cell:'
        write(*,*)ixmin,iymin,izmin
        write(*,*)'To cell:'
        write(*,*)ixmax,iymax,izmax

! set the slot volume to free space
    
        do ix=ixmin,ixmax
          do iy=iymin,iymax
            do iz=izmin,izmax
              material_mesh(centre,ix,iy,iz)=0
            end do
          end do
        end do
        
        
      end do  ! next slot
    
    else
   
      write(*,*)'ERROR: unknown additional component type:',trim(component_type)
    
    end if
  
  end do ! next component

  write(*,*)'FINISHED specify_additional_components'

RETURN  
  
END SUBROUTINE specify_additional_components
