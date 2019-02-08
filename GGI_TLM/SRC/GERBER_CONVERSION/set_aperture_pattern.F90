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
!     SUBROUTINE set_aperture_pattern
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
SUBROUTINE set_aperture_pattern(rloop)

USE gerber

IMPLICIT NONE

! local_variables

integer :: rloop

integer iax,iay
real*8  :: lx,ly,rx,ry
real*8  :: ax,ay,ar
integer :: ns
real*8  :: theta,theta_0,alpha

real*8 :: xp,yp,xn,yn,pline,ptest
integer :: i

real*8,parameter :: pi=3.1415926535d0

! set the aperture pattern

! START

! only set the aperture on the second time through
  if (rloop.eq.1) RETURN

  if (allocated( ap )) Deallocate( ap )

  anx=NINT( 0.5d0*Asize(aperture)/dl )+1
  any=NINT( 0.5d0*Asize(aperture)/dl )+1
 
  write(*,*)'Aperture size is:',Asize(aperture)
  write(*,*)'dl=',dl
  write(*,*)'anx=',anx,' any=',any
 
  ALLOCATE( ap(-anx:anx,-any:any) )
  ap(-anx:anx,-any:any)=0
 
! loop over the aperture points        
  do iax=-anx,anx
    do iay=-any,any  
! get the coordinates of this point on the aperture
  
      ax=iax*dl
      ay=iay*dl
  
      if (Atype(aperture).EQ.'C') then
! Circular aperture with hole
    
        ar=sqrt(ax*ax+ay*ay)
        if ( (ar.LE.Aparams(aperture,1)/2d0).AND.(ar.GE.Aparams(aperture,2)/2d0) ) then
          ap(iax,iay)=1
        else if ( ar.LT.Aparams(aperture,2)/2d0 ) then
          ap(iax,iay)=0
        end if
      
      else if (Atype(aperture).EQ.'R') then
! Rectangular aperture with hole
      
        ar=sqrt(ax*ax+ay*ay)
        if ( (abs(ax).LE.Aparams(aperture,1)/2d0).AND.   &
             (abs(ay).LE.Aparams(aperture,2)/2d0).AND.   &
             (ar.GE.Aparams(aperture,3)/2d0) ) then
          ap(iax,iay)=1
        else if ( ar.LT.Aparams(aperture,3)/2d0 ) then
          ap(iax,iay)=0
        end if

      else if (Atype(aperture).EQ.'O') then
! Obround aperture with hole

        lx=Aparams(aperture,1)/2d0
        ly=Aparams(aperture,2)/2d0
        
! set rectangle first
        if (lx.gt.ly) then
          rx=lx-ly
          ry=ly
          if ( (abs(ax).LE.rx).AND.(abs(ay).LE.ry) ) ap(iax,iay)=1
! LH semi circle
          ar=sqrt((ax-rx)**2+ay*ay)
          if ( ar.LT.ly ) ap(iax,iay)=1
          
! RH semicircle
          ar=sqrt((ax+rx)**2+ay*ay)
          if ( ar.LT.ly ) ap(iax,iay)=1
          
        else
          rx=lx
          ry=ly-lx
          if ( (abs(ax).LE.rx).AND.(abs(ay).LE.ry) ) ap(iax,iay)=1
! top semi circle
          ar=sqrt(ax*ax+(ay-ry)**2)
          if ( ar.LT.lx ) ap(iax,iay)=1
          
! bottom semicircle
          ar=sqrt(ax*ax+(ay+ry)**2)
          if ( ar.LT.lx ) ap(iax,iay)=1
       
        end if
  
! hole
        ar=sqrt(ax*ax+ay*ay)
        if ( ar.LT.Aparams(aperture,3)/2d0 ) then
          ap(iax,iay)=0
        end if
      
  
      else if (Atype(aperture).EQ.'P') then
               
! enclosing circle
        ar=sqrt(ax*ax+ay*ay)
        if ( ar.LT.Aparams(aperture,1)/2d0 ) then
          ap(iax,iay)=1
        end if
        
        ns=NINT(Aparams(aperture,2))           ! number of sides
        
        if (ns.LT.3) then
          write(*,*)'Number of sides in polygon is less than 3'
          write(*,*)'ns=',ns,Aparams(aperture,2)
          STOP
        else if (ns.GT.12) then
          write(*,*)'Number of sides in polygon is greater than 12'
          write(*,*)'ns=',ns,Aparams(aperture,2)
          STOP
        end if
        
        theta_0=Aparams(aperture,3)*pi/180d0   ! offset angle
        
! loop over sides
        do i=1,ns
        
          theta=theta_0+(dble(i)-0.5d0)*2d0*pi/dble(ns) ! angle to the midpoint of the edge
          xn=cos(theta)  ! normal to the edge
          yn=sin(theta)
          
          alpha=theta_0+dble(i-1)*2d0*pi/dble(ns)  ! angle to a point on the edge (first vertex)
          
          xp=(Aparams(aperture,1)/2d0)*cos(alpha)  ! get a point on the edge (first vertex)
          yp=(Aparams(aperture,1)/2d0)*sin(alpha)
          
          pline=xp*xn+yp*yn   ! point on edge dotted with normal direction
          
          ptest=ax*xn+ay*yn   ! point dotted with normal direction
          if (ptest.GT.pline) ap(iax,iay)=0
        
        end do
        
! hole
        ar=sqrt(ax*ax+ay*ay)
        if ( ar.LT.Aparams(aperture,4)/2d0 ) then
          ap(iax,iay)=0
        end if

      end if
  
    end do
  end do        
 
  RETURN
          
END SUBROUTINE set_aperture_pattern
