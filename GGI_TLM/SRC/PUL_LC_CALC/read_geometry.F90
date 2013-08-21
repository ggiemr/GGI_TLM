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
  SUBROUTINE read_geometry()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  

  integer wire
  integer pul_nterms_in
  integer row
  character ch
  character*80 ipfile
    
! START

  write(*,*)'Enter filename for geometry'
  read(*,'(A80)')ipfile
  open(unit=8,file=ipfile)
  
  read(8,*)pul_nwires_in
  read(8,*)pul_dl
  
  read(8,'(A1)')ch
  if( (ch.eq.'y').OR.(ch.eq.'Y') ) then
    pul_include_dielectric=.TRUE.
  else 
    pul_include_dielectric=.FALSE.
  end if
  write(*,*)'Include dielectric:',pul_include_dielectric
  
  read(8,'(A1)')ch
  if( (ch.eq.'y').OR.(ch.eq.'Y') ) then
    pul_include_TLM_return=.TRUE.
  else 
    pul_include_TLM_return=.FALSE.
  end if
  write(*,*)'Include TLM return:',pul_include_TLM_return
  
  if (pul_include_TLM_return) then
    pul_nwires=pul_nwires_in+1
  else
    pul_nwires=pul_nwires_in
  end if
  
  read(8,*)pul_return_conductor
  
  read(8,*)pul_op_flag
  read(8,*)pul_nterms_in

! allocate memory for wire specification

  write(*,*)'Allocating wire_spec'
  allocate ( pul_wire_spec(1:pul_nwires) )
  
! allocate individual wire data

  read(8,*)

  do wire=1,pul_nwires_in
  
    pul_wire_spec(wire)%nterms=pul_nterms_in
    pul_wire_spec(wire)%npoints=1+2*pul_wire_spec(wire)%nterms
    
    read(8,*)row,pul_wire_spec(wire)%xc  &
                ,pul_wire_spec(wire)%yc  &
                ,pul_wire_spec(wire)%rw  &
                ,pul_wire_spec(wire)%ri  &
                ,pul_wire_spec(wire)%epsr 
        
    if (row.ne.wire) then
      write(*,*)'Wire data should be provided in order'
      write(*,*)'Expecting data for wire number ',wire
      stop
    end if 	
	
  end do ! next wire
  
  if (pul_include_TLM_return) then
    wire=pul_nwires
    pul_wire_spec(wire)%xc=0d0
    pul_wire_spec(wire)%yc=0d0
    pul_wire_spec(wire)%rw=pul_dl*1.08d0/2d0 
    pul_wire_spec(wire)%ri=pul_wire_spec(wire)%rw
    pul_wire_spec(wire)%epsr=1.0
    pul_wire_spec(wire)%nterms=pul_nterms_in
    pul_wire_spec(wire)%npoints=1+2*pul_wire_spec(wire)%nterms
  end if
  
  read(8,*)pul_xmin,pul_xmax,pul_nxp
  read(8,*)pul_ymin,pul_ymax,pul_nyp
  
  close(unit=8)  

  return
  
  end SUBROUTINE read_geometry
  
