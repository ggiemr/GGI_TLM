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
!     started 26/3/19 CJS   add dielectric volumes 
!     
!

SUBROUTINE specify_additional_components

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i

character(LEN=256) :: component_type
character(LEN=256) :: material_name

real :: hxmin,hxmax,hymin,hymax,hzmin,hzmax
integer :: nslots
real*8  :: ws,wr,ds

character(LEN=1) :: slot_direction
character(LEN=1) :: slot_normal
character(LEN=4) :: slot_face

integer :: normx,normy,normz

integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
real*8  :: lx,ly,lz,l
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
    
      lx=hxmax-hxmin
      ly=hymax-hymin
      lz=hzmax-hzmin
    
      read(10,*)nslots
      
      slot_xmin=hxmin
      slot_xmax=hxmax
      slot_ymin=hymin
      slot_ymax=hymax
      slot_zmin=hzmin
      slot_zmax=hzmax
      normx=0
      normy=0
      normz=0
      
      read(10,'(A)')slot_direction
      
      if (slot_direction.EQ.'x') then
        normx=1
      else if (slot_direction.EQ.'y') then
        normy=1
      else if (slot_direction.EQ.'z') then
        normz=1
      else
        write(*,*)'Heatsink slot direction should be x, y or z'
        STOP 1
      end if
      
      read(10,'(A4)')slot_face
      
      read(10,*)ws
      
      read(10,*)ds
            
      if      (slot_face.EQ.'xmin') then
        slot_xmax=hxmin+ds                   
        normx=1
      else if (slot_face.EQ.'xmax') then
        slot_xmin=hxmax-ds                
        normx=1
      else if (slot_face.EQ.'ymin') then
        slot_ymax=hymin+ds                 
        normy=1
      else if (slot_face.EQ.'ymax') then
        slot_ymin=hymax-ds                
        normy=1
      else if (slot_face.EQ.'zmin') then
        slot_zmax=hzmin+ds                
        normz=1
      else if (slot_face.EQ.'zmax') then
        slot_zmin=hzmax-ds                        
        normz=1
      else
      
        write(*,*)'ERROR: unknown slot_face:',slot_face
        STOP 1
        
      end if

! work out the slot_normal direction

      if (normx.Eq.0) then
        slot_normal='x'
        l=lx
      else if (normy.Eq.0) then
        slot_normal='y'
        l=ly
      else if (normz.Eq.0) then
        slot_normal='z'
        l=lz     
      else
        write(*,*)'ERROR specifying the slot normal direction'
        STOP 1
      end if  
    
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

      wr=(l-nslots*ws)/dble(nslots+1)
      if (wr.LT.0d0) then
        write(*,*)'ERROR: width of ribs in heatsink is less than zero'
        write(*,*)'wr=',wr
        STOP 1
      end if
      
      write(*,*)'slot_direction:',slot_direction
      write(*,*)'slot_face:',slot_face
      write(*,*)'slot_normal:',slot_normal
      
      write(*,*)'slot_xmin=',slot_xmin,' slot_xmax=',slot_xmax
      write(*,*)'slot_ymin=',slot_ymin,' slot_ymax=',slot_ymax
      write(*,*)'slot_zmin=',slot_zmin,' slot_zmax=',slot_zmax
      
      write(*,*)'Heatsink slot width:',ws
      write(*,*)'Heatsink rib width :',wr
      write(*,*)'Heatsink depth     :',ds

! assume slots are in the zmin face and are in the y direction

      do slot=1,nslots
      
        if (slot_normal.Eq.'x') then
          slot_xmin=hxmin+wr+(slot-1)*(wr+ws)
          slot_xmax=slot_xmin+ws
        else if (slot_normal.Eq.'y') then
          slot_ymin=hymin+wr+(slot-1)*(wr+ws)
          slot_ymax=slot_ymin+ws
        else if (slot_normal.Eq.'z') then
          slot_zmin=hzmin+wr+(slot-1)*(wr+ws)
          slot_zmax=slot_zmin+ws        
        end if
        
        write(*,*)'Slot number:',slot
        write(*,*)'slot_xmin=',slot_xmin,' slot_xmax=',slot_xmax
        write(*,*)'slot_ymin=',slot_ymin,' slot_ymax=',slot_ymax
        write(*,*)'slot_zmin=',slot_zmin,' slot_zmax=',slot_zmax
        
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
    
    else if( index(component_type,'dielectric').NE.0 ) then
    
      read(10,*)hxmin,hymin,hzmin,hxmax,hymax,hzmax
    
      read(10,'(A)')material_name

      n_volumes=n_volumes+1
    
      if(n_volumes.GT.max_volumes) then
    
        write(*,*)'ERROR in GGI_TLM_create_PCB_simulation_model: maximum number of volumes exceeded'
        write(*,*)'Maximum number of volumes is set to ',max_volumes
        write(*,*)'in /GGI_TLM/SRC/TLM_MODULES/PCB_simulation_setup_modules.F90'
        STOP 1
    
      end if

      volume_type(n_volumes)=volume_type_rectangular_block2
    
      volume_parameters(n_volumes,1)=hxmin
      volume_parameters(n_volumes,2)=hymin
      volume_parameters(n_volumes,3)=hzmin
      volume_parameters(n_volumes,4)=hxmax
      volume_parameters(n_volumes,5)=hymax
      volume_parameters(n_volumes,6)=hzmax

      n_volume_materials=n_volume_materials+1
      volume_material_type(n_volume_materials)=volume_material_type_DISPERSIVE
      volume_material_name(n_volume_materials)=trim(material_name)
      volume_material_to_volume_list(n_volume_materials)=n_volumes         
    
    else
   
      write(*,*)'ERROR: unknown additional component type:',trim(component_type)
    
    end if
  
  end do ! next component

  write(*,*)'FINISHED specify_additional_components'

RETURN  
  
END SUBROUTINE specify_additional_components
