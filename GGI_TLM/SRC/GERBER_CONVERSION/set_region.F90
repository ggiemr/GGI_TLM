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
!     SUBROUTINE set_region
!     SUBROUTINE fill_4
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 10/12/19 CJS
!     
!
SUBROUTINE set_region(r,p,nx,ny,polarity)

USE gerber_parameters

IMPLICIT NONE

integer :: nx,ny
integer :: r(1:nx,1:ny)
integer :: p(1:nx,1:ny)
integer :: polarity

integer :: ix,iy

! START

  write(*,*)'CALLED set_region'

! step around the outer edge looking for zeros. 
! If found, call the fill routine
  write(*,*)'Fill exterior region'
  do ix=1,nx
    if (r(ix,1).EQ.0)  CALL fill_4(nx,ny,r,ix,1,0,-1)
    if (r(ix,ny).EQ.0) CALL fill_4(nx,ny,r,ix,nx,0,-1)
  end do
  do iy=1,ny
    if (r(1,iy).EQ.0)  CALL fill_4(nx,ny,r,1,iy,0,-1)
    if (r(nx,iy).EQ.0) CALL fill_4(nx,ny,r,nx,iy,0,-1)
  end do
  
! we have now filled the exterior region
  write(*,*)'Fill defined regions'

! we now loop over the area looking for unfilled regions and fill them
  do ix=1,nx
    do iy=1,ny
      if (r(ix,iy).EQ.0)  then
        CALL fill_4(nx,ny,r,ix,iy,0,1)
      end if
    end do 
  end do
  
! transfer the region fill array to the pixel array
  write(*,*)'transfer the region fill array to the pixel array'  
  do ix=1,nx
    do iy=1,ny
        
      if (polarity.EQ.dark) then
        if (r(ix,iy).EQ.1) then
          p(ix,iy)=1
        end if
      else if (polarity.EQ.clear) then
        if (r(ix,iy).EQ.1) then
          p(ix,iy)=0
        end if        
      end if

    end do 
  end do
  
! reset the region fill array

  r(1:nx,1:ny)=0

  write(*,*)'Done: set_region'

END SUBROUTINE set_region
!
! ___________________________________________________________________
!
!
SUBROUTINE fill_4(nx,ny,p,px,py,n0,n1)

! flood fill cells with value n0 and replace with n1

integer :: nx,ny
integer :: p(nx,ny)
integer :: px,py,n0,n1

integer :: len_q
integer :: q(1:nx*ny,2)
integer :: ix,iy

! START

! return if the point is out of range
if ( (px.LT.1).OR.(px.GT.nx).OR.(py.LT.1).OR.(py.GT.ny) ) RETURN

! return if the point is not the value for replacement
if (p(px,py).NE.n0) RETURN

len_Q=1
q(len_Q,1)=px
q(len_Q,2)=py

do while (len_q.GE.1)

! change the value of the current point
  ix=q(len_Q,1)
  iy=q(len_Q,2)
  p(ix,iy)=n1
  
! remove this point from the list
  len_Q=len_Q-1
  
! examine the surrounding points and add to the queue 
  if ((ix-1).GE.1) then
    if (p(ix-1,iy).EQ.n0) then
      len_Q=len_Q+1
      if(len_Q.GT.nx*ny) GOTO 9000
      q(len_Q,1)=ix-1
      q(len_Q,2)=iy
    end if
  end if
  
  if ((ix+1).LE.nx) then
    if (p(ix+1,iy).EQ.n0) then
      len_Q=len_Q+1
      if(len_Q.GT.nx*ny) GOTO 9000
      q(len_Q,1)=ix+1
      q(len_Q,2)=iy
    end if
  end if
  
  if ((iy-1).GE.1) then
    if (p(ix,iy-1).EQ.n0) then
      len_Q=len_Q+1
      if(len_Q.GT.nx*ny) GOTO 9000
      q(len_Q,1)=ix
      q(len_Q,2)=iy-1
    end if
  end if
  
  if ((iy+1).LE.ny) then
    if (p(ix,iy+1).EQ.n0) then
      len_Q=len_Q+1
      if(len_Q.GT.nx*ny) GOTO 9000
      q(len_Q,1)=ix
      q(len_Q,2)=iy+1
    end if
  end if
 
end do

RETURN

9000 write(*,*)'ERROR: queue size exceeded in fill_4'
STOP

END 

