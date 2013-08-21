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
  SUBROUTINE calc_dielectric_points()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  

  real*8 x,y,r,t
  real*8 xp,yp,rp,tp,sp

  integer wire
  integer tloop
  
! START
  
  if (pul_include_dielectric) then
!    write(*,*)'Calculate dielectric points'

! allocate dielecric points and set positions

    if (pul_op_flag)    write(10,*)'Dielectric points'
    if (pul_op_flag)    write(10,*)  &
'wire     xc         yc       tw         xp       yp         nx        ny '

    do wire=1,pul_nwires
    
      if (allocated( pul_wire_spec(wire)%xdp ) ) deallocate ( pul_wire_spec(wire)%xdp )
      if (allocated( pul_wire_spec(wire)%ydp ) ) deallocate ( pul_wire_spec(wire)%ydp )
      if (allocated( pul_wire_spec(wire)%rdp ) ) deallocate ( pul_wire_spec(wire)%rdp )
      if (allocated( pul_wire_spec(wire)%tdp ) ) deallocate ( pul_wire_spec(wire)%tdp )
      
      if (allocated( pul_wire_spec(wire)%nxdp ) ) deallocate ( pul_wire_spec(wire)%nxdp )
      if (allocated( pul_wire_spec(wire)%nydp ) ) deallocate ( pul_wire_spec(wire)%nydp )
     
      if (allocated( pul_wire_spec(wire)%a2 ) ) deallocate ( pul_wire_spec(wire)%a2 )
      if (allocated( pul_wire_spec(wire)%b2 ) ) deallocate ( pul_wire_spec(wire)%b2 )
       
      allocate ( pul_wire_spec(wire)%xdp(1:pul_wire_spec(wire)%npoints) )
      allocate ( pul_wire_spec(wire)%ydp(1:pul_wire_spec(wire)%npoints) )
      allocate ( pul_wire_spec(wire)%rdp(1:pul_wire_spec(wire)%npoints) )
      allocate ( pul_wire_spec(wire)%tdp(1:pul_wire_spec(wire)%npoints) )
      
      allocate ( pul_wire_spec(wire)%nxdp(1:pul_wire_spec(wire)%npoints) )
      allocate ( pul_wire_spec(wire)%nydp(1:pul_wire_spec(wire)%npoints) )
    
      allocate ( pul_wire_spec(wire)%a2(1:pul_wire_spec(wire)%npoints) )
      allocate ( pul_wire_spec(wire)%b2(1:pul_wire_spec(wire)%npoints) )
! wire center coordinates    

      x=pul_wire_spec(wire)%xc
      y=pul_wire_spec(wire)%yc
      r=pul_wire_spec(wire)%ri
    
      if (r.eq.0d0) then
        write(*,*)'Error in calc_dielectric_points'
        write(*,*)'Dielectric radius is zero'
        STOP
      end if

! loop over theta    

      do tloop=1,pul_wire_spec(wire)%npoints
    
        t=2d0*pi*dble(tloop-1)/dble(pul_wire_spec(wire)%npoints)
        xp=x+r*cos(t)
        yp=y+r*sin(t)
        rp=sqrt(xp*xp+yp*yp)
        tp=atan2(yp,xp)
        pul_wire_spec(wire)%xdp(tloop)=xp
        pul_wire_spec(wire)%ydp(tloop)=yp
        pul_wire_spec(wire)%rdp(tloop)=rp
        pul_wire_spec(wire)%tdp(tloop)=tp

! normal to dielectric surface at point	
        pul_wire_spec(wire)%nxdp(tloop)=cos(t)
        pul_wire_spec(wire)%nydp(tloop)=sin(t)
      
!        write(*,8000)wire,x,y,r,xp,yp,rp,tp
8000    format(I4,7F10.6)

        if (pul_op_flag) write(10,8000)wire,x,y,r,xp,yp,pul_wire_spec(wire)%nxdp(tloop),pul_wire_spec(wire)%nydp(tloop)
      
      end do
    
    end do ! next wire

  end if  ! include dielectric
  
  return
  
  end SUBROUTINE calc_dielectric_points
  
