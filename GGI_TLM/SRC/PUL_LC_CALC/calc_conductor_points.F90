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
  SUBROUTINE calc_conductor_points()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  

  real*8 x,y,r,t
  real*8 xp,yp,rp,tp,sp

  integer wire
  integer tloop
  
! START
  
!  write(*,*)'CALLED: Calculate conductor points'

! allocate wire points and set positions
  
!    if (pul_op_flag)    write(10,*)'Conductor points'
!  
!    if (pul_op_flag)    write(10,*)  &
!'wire     xc         yc       tw         xp       yp         rp        tp '
    
  do wire=1,pul_nwires
      
    if (allocated( pul_wire_spec(wire)%rp ) ) then
      deallocate ( pul_wire_spec(wire)%rp )
    end if
      
    if (allocated( pul_wire_spec(wire)%tp ) ) then
      deallocate ( pul_wire_spec(wire)%tp )
    end if
    
    if (allocated( pul_wire_spec(wire)%a ) )then
      deallocate ( pul_wire_spec(wire)%a )
    end if
       
    if (allocated( pul_wire_spec(wire)%b ) )then
      deallocate ( pul_wire_spec(wire)%b )
    end if
    
    if (allocated( pul_wire_spec(wire)%yp ) ) then
      deallocate ( pul_wire_spec(wire)%yp )
    end if
    
    if (allocated( pul_wire_spec(wire)%xp ) ) then
      deallocate ( pul_wire_spec(wire)%xp )
    end if
  
    allocate ( pul_wire_spec(wire)%xp(1:pul_wire_spec(wire)%npoints) )
    allocate ( pul_wire_spec(wire)%yp(1:pul_wire_spec(wire)%npoints) )
    allocate ( pul_wire_spec(wire)%rp(1:pul_wire_spec(wire)%npoints) )
    allocate ( pul_wire_spec(wire)%tp(1:pul_wire_spec(wire)%npoints) )
    
    allocate ( pul_wire_spec(wire)%a(1:pul_wire_spec(wire)%npoints) )
    allocate ( pul_wire_spec(wire)%b(1:pul_wire_spec(wire)%npoints) )
    
! wire center coordinates    

    x=pul_wire_spec(wire)%xc
    y=pul_wire_spec(wire)%yc
    r=pul_wire_spec(wire)%rw
    
    if (r.eq.0d0) then
      write(*,*)'Error in calc_conductor_points'
      write(*,*)'Conductor radius is zero'
      STOP
    end if
    
! loop over theta    
    do tloop=1,pul_wire_spec(wire)%npoints
    
      t=2d0*pi*dble(tloop-1)/dble(pul_wire_spec(wire)%npoints)
      xp=x+r*cos(t)
      yp=y+r*sin(t)
      rp=sqrt(xp*xp+yp*yp)
      tp=atan2(yp,xp)
      pul_wire_spec(wire)%xp(tloop)=xp
      pul_wire_spec(wire)%yp(tloop)=yp
      pul_wire_spec(wire)%rp(tloop)=rp
      pul_wire_spec(wire)%tp(tloop)=tp
      
      if (pul_op_flag) write(10,8000)wire,x,y,r,xp,yp,rp,tp
8000  format(I4,7F10.6)
      
    end do
    
  end do ! next wire
  
!  write(*,*)'FINISHED: Calculate conductor points'
  
  return
  
  end SUBROUTINE calc_conductor_points
  
