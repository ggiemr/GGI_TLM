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
!     SUBROUTINE sweep_line
!     SUBROUTINE get_r_theta
!     SUBROUTINE get_centre_clockwise
!     SUBROUTINE get_centre_counterclockwise
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
SUBROUTINE sweep_line(p,nx,ny,ap,anx,any,x,y,xm,ym,dx,dy,xmin,ymin,dl,interpolation_mode,quadrant_mode,polarity,region)

USE gerber_parameters
USE constants

IMPLICIT NONE

integer :: nx,ny,anx,any
integer :: p(1:nx,1:ny)
integer :: ap(-anx:anx,-any:any)
real*8  :: x,y,xm,ym,dx,dy,xmin,ymin,dl
integer :: interpolation_mode,quadrant_mode
integer :: polarity
logical :: region

integer :: ipx,ipy,iax,iay,ix,iy
real*8  :: sx,sy,dsx,dsy
real*8  :: length
integer :: n,i

real*8  :: xi,yj
real*8  :: t1,t2,dt,r1,r2,dr,r,t

! START
  if ( (interpolation_mode.NE.linear).AND.     &
       (interpolation_mode.NE.clockwise).AND.     &
       (interpolation_mode.NE.counterclockwise) ) then
  
    write(*,*)'Unknown interpolation mode:',interpolation_mode
    STOP
    
  end if

  if (interpolation_mode.EQ.linear) then
  
    write(*,*)'Linear'
    
! set the spatial step for the linear sweep
    
    length=sqrt( (xm-x)**2+(ym-y)**2 )
    if (length.GT.0d0) then
      n=NINT( length/(0.5d0*dl) )
    else
      n=1
    end if
    n=max(n,1)
    
! If we are going to define a region then make sure we have a small enough spatial step to ensure no holes in the line
    if (region) then
      n=n*2.0
    end if
    
    dsx=(xm-x)/n
    dsy=(ym-y)/n
            
! flash the aperture at points along the line  

    do i=0,n

      sx=x+dsx*i
      sy=y+dsy*i      
      if (region) then
        ipx=NINT((sx-xmin)/dl)
        ipy=NINT((sy-ymin)/dl)
        if ((ipx.GT.0).AND.(ipx.LE.nx).AND.(ipy.GT.0).AND.(ipy.LE.ny)) then
          p(ipx,ipy)=1
        end if
      else
        CALL flash(p,nx,ny,ap,anx,any,sx,sy,xmin,ymin,dl,polarity)
      end if
         
    end do
    
    RETURN

  end if

! must be a circular arc...
  
  if (interpolation_mode.EQ.clockwise) then
  
    if (quadrant_mode.EQ.single) then
      write(*,*)'Clockwise, single'
  
      CALL get_centre_clockwise(xm,ym,x,y,dx,dy,xi,yj,r1,r2,t1,t2)
      
    else if (quadrant_mode.EQ.multi) then
      write(*,*)'Clockwise, multi'
  
      xi=xm+dx
      yj=ym+dy
      
      CALL get_r_theta(xm,ym,x,y,xi,yj,r1,r2,t1,t2)

! ensure that t1 > t2 i.e. we go clockwise from t1 to t2
      if (t1.LE.t2) t1=t1+2d0*pi 
      
    end if   

  else if (interpolation_mode.EQ.counterclockwise) then
  
    if (quadrant_mode.EQ.single) then
      write(*,*)'Anticlockwise, single'
    
      CALL get_centre_counterclockwise(xm,ym,x,y,dx,dy,xi,yj,r1,r2,t1,t2)
      
    else if (quadrant_mode.EQ.multi) then
      write(*,*)'Anticlockwise, multi'
  
      xi=xm+dx
      yj=ym+dy
      
      CALL get_r_theta(xm,ym,x,y,xi,yj,r1,r2,t1,t2)
      
! ensure that t2 > t1 i.e. we go counterclockwise from t1 to t2
      if (t2.LE.t1) t2=t2+2d0*pi 
      
    end if   
    
  end if
  
! Sweep along the arc flashing the aperture

! set the spatial step for the arc sweep
    
  length=((r2+r1)/2d0)*abs(t2-t1)    
  if (length.GT.0d0) then
    n=NINT( length/(0.5d0*dl) )
  else
    n=1
  end if
  n=max(n,1)
    
! If we are going to define a region then make sure we have a small enough spatial step to ensure no holes in the line
  if (region) then
    n=n*2.0
  end if
    
  dt=(t2-t1)/n
  dr=(r2-r1)/n
  
  write(*,*)'length=',length,' n=',n
  write(*,*)'r1=',r1,' r2=',r2,' dr=',dr
  write(*,*)'theta1=',t1*180d0/pi,' theta2=',t2*180d0/pi,' dt=',dt
  write(*,*)'centre=',xi,yj
            
! flash the aperture at points along the arc 

  do i=0,n

    t=t1+dt*i
    r=r1+dr*i
    sx=r*cos(t)+xi
    sy=r*sin(t)+yj  
    if (region) then
      ipx=NINT((sx-xmin)/dl)
      ipy=NINT((sy-ymin)/dl)
      if ((ipx.GT.0).AND.(ipx.LE.nx).AND.(ipy.GT.0).AND.(ipy.LE.ny)) then
        p(ipx,ipy)=1
      end if
    else
      CALL flash(p,nx,ny,ap,anx,any,sx,sy,xmin,ymin,dl,polarity)
    end if
       
  end do

  
  write(*,*)'Done: sweep_line'

END SUBROUTINE sweep_line
!
! ___________________________________________________________
!
!
SUBROUTINE get_r_theta(xm,ym,x,y,xi,yj,r1,r2,t1,t2)

USE constants

IMPLICIT NONE

real*8 xm,ym,x,y,xi,yj,r1,r2,t1,t2

! START
 
 r1=sqrt((xi-xm)**2+(yj-ym)**2)
 r2=sqrt((xi-x)**2+(yj-y)**2)
 if ( (ym-yj.Eq.0d0).AND.(xm-xi.Eq.0d0) ) then
   write(*,*)'0/0 Error in atan2'
 end if
 t1=atan2(ym-yj,xm-xi)
 if (t1.LT.0d0) t1=t1+2d0*pi  ! ensure positive angles
 
 if ( (y-yj.Eq.0d0).AND.(x-xi.Eq.0d0) ) then
   write(*,*)'0/0 Error in atan2'
 end if
 t2=atan2(y-yj,x-xi)
 if (t2.LT.0d0) t2=t2+2d0*pi  ! ensure positive angles


END SUBROUTINE get_r_theta
!
! ___________________________________________________________
!
!
SUBROUTINE get_centre_clockwise(xm,ym,x,y,dx,dy,xi,yj,r1,r2,t1,t2)

USE constants

IMPLICIT NONE

real*8 xm,ym,x,y,dx,dy,xi,yj,r1,r2,t1,t2

real*8 :: sign(4,2)
real*8 :: dr(4) ! radius deviation for each centre
real*8 :: min_dr
real*8 :: dt(4) ! theta deviation for each centre
integer :: i,min_i

! START

sign(1,1)=1
sign(1,2)=1

sign(2,1)=1
sign(2,2)=-1

sign(3,1)=-1
sign(3,2)=1

sign(4,1)=-1
sign(4,2)=-1

min_dr=1d30
min_i=0

do i=1,4

  xi=xm+sign(i,1)*dx
  yj=ym+sign(i,2)*dy
  
  CALL get_r_theta(xm,ym,x,y,xi,yj,r1,r2,t1,t2)
  
! ensure that t1 > t2 i.e. we go clockwise from t1 to t2
  do while (t1.LE.t2) 
    t1=t1+2d0*pi
  end do
  
  dr(i)=r2-r1
  dt(i)=t2-t1
  
  if (dt(i).GE.-(pi/2d0)*1.0001d0) then
    if (dr(i).LT.min_dr) then
      min_dr=min(min_dr,dr(i))
      min_i=i
    end if
  end if
  
end do

if (min_i.EQ.0) then
  write(*,*)'Error in get_centre_clockwise: cannot find a viable centre'
  write(*,*)'xm,ym:',xm,ym
  write(*,*)'x,y  :',x,y
  write(*,*)'dx,dy:',dx,dy
  do i=1,4
    write(*,'(4ES16.6)')xi,yj,dr(i),dt(i)
  end do
  STOP
end if


xi=xm+sign(min_i,1)*dx
yj=ym+sign(min_i,2)*dy

CALL get_r_theta(xm,ym,x,y,xi,yj,r1,r2,t1,t2)
  
! ensure that t1 > t2 i.e. we go clockwise from t1 to t2
do while (t1.LE.t2) 
  t1=t1+2d0*pi
end do

RETURN

END SUBROUTINE get_centre_clockwise
!
! ___________________________________________________________
!
!
SUBROUTINE get_centre_counterclockwise(xm,ym,x,y,dx,dy,xi,yj,r1,r2,t1,t2)

USE constants

IMPLICIT NONE

real*8 xm,ym,x,y,dx,dy,xi,yj,r1,r2,t1,t2

real*8 :: sign(4,2)
real*8 :: dr(4) ! radius deviation for each centre
real*8 :: min_dr
real*8 :: dt(4) ! theta deviation for each centre
integer :: i,min_i

! START

sign(1,1)=1
sign(1,2)=1

sign(2,1)=1
sign(2,2)=-1

sign(3,1)=-1
sign(3,2)=1

sign(4,1)=-1
sign(4,2)=-1

min_dr=1d30
min_i=0

do i=1,4

  xi=xm+sign(i,1)*dx
  yj=ym+sign(i,2)*dy
  
  CALL get_r_theta(xm,ym,x,y,xi,yj,r1,r2,t1,t2)
  
! ensure that t1<t2 i.e. we go counterclockwise from t1 to t2
  do while (t2.LE.t1) 
    t2=t2+2d0*pi
  end do
  
  write(*,*)'t2>t1',t2*180d0/pi,t1*180d0/pi
  
  dr(i)=r2-r1
  dt(i)=t2-t1
  
  if (dt(i).LE.(pi/2d0)*1.0001d0) then
    if (dr(i).LT.min_dr) then
      min_dr=min(min_dr,dr(i))
      min_i=i
    end if
  end if
  
end do

if (min_i.EQ.0) then
  write(*,*)'Error in get_centre_clockwise: cannot find a viable centre'
  write(*,*)'xm,ym:',xm,ym
  write(*,*)'x,y  :',x,y
  write(*,*)'dx,dy:',dx,dy
  do i=1,4
    write(*,'(4ES16.6)')xi,yj,dr(i),dt(i)
  end do
  STOP
end if


xi=xm+sign(min_i,1)*dx
yj=ym+sign(min_i,2)*dy

CALL get_r_theta(xm,ym,x,y,xi,yj,r1,r2,t1,t2)
  
! ensure that t1<t2 i.e. we go counterclockwise from t1 to t2
do while (t2.LE.t1) 
  t2=t2+2d0*pi
end do

RETURN

END SUBROUTINE get_centre_counterclockwise
