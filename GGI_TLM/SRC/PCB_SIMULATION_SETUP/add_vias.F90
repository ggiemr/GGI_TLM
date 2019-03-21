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
!     SUBROUTINE add_vias
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

SUBROUTINE add_vias

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i
real*8  :: vxmin,vymin,vzmin        ! minimum coordinates of via
real*8  :: vxmax,vymax,vzmax        ! maximum coordinates of via

! START

  write(*,*)'CALLED: add_vias'

  read(10,*)n_vias
  
  do i=1,n_vias
   
    read(10,*)vxmin,vymin,vzmin,vzmax 
    
    vxmax=vxmin
    vymax=vymin
    
    write(*,*)'min coordinate of via:'
    write(*,*)vxmin,vymin,vzmin
    write(*,*)'max coordinate of via:'
    write(*,*)vxmax,vymax,vzmax

! Set the cells in the mesh. Note that we offset the z coordinates slightly in case the z value is on a cell face 
! this leads to ambiguity as to which cell the end point is situated in

    CALL set_terminal_connection_cells(vxmin,vymin,vzmin+dl/20,vxmax,vymax,vzmax-dl/20)

  end do

  write(*,*)'FINISHED: add_vias'

RETURN  
  
END SUBROUTINE add_vias
