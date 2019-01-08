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
!     SUBROUTINE flash
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
SUBROUTINE flash(p,nx,ny,ap,anx,any,x,y,xmin,ymin,dl,polarity)

USE gerber_parameters

IMPLICIT NONE

integer :: nx,ny,anx,any
integer :: p(1:nx,1:ny)
integer :: ap(-anx:anx,-any:any)
real*8  :: x,y,xmin,ymin,dl
integer :: polarity

integer :: ipx,ipy,iax,iay,ix,iy

! START

! get the centre point in the pixel array 

  ipx=NINT((x-xmin)/dl)
  ipy=NINT((y-ymin)/dl)
  
! loop over the aperture points        
  do iax=-anx,anx
    do iay=-any,any  
  
      ix=ipx+iax
      iy=ipy+iay
    
      if ((ix.GT.0).AND.(ix.LE.nx).AND.(iy.GT.0).AND.(iy.LE.ny)) then
        
        if (polarity.EQ.dark) then
          if (ap(iax,iay).EQ.1) then
            p(ix,iy)=1
          end if
        else if (polarity.EQ.clear) then
          if (ap(iax,iay).EQ.1) then
            p(ix,iy)=0
          end if        
        end if
        
      end if
  
    end do
  end do        

END SUBROUTINE flash
