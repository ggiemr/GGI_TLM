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
  SUBROUTINE write_geometry()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  

  real*8 x,y,r,t
  real*8 xp,yp,rp,tp,sp

  integer wire
  integer tloop
  
! START
  
! open files for output test data

  if (pul_op_flag) then
  
    write(*,*)'Write geometry'
  
    open(unit=10,file='pul_LC_points.dat')  
    
    open(unit=11,file='pul_LC_cell.dat')  
    write(11,*)-pul_dl/2d0,-pul_dl/2d0
    write(11,*) pul_dl/2d0,-pul_dl/2d0
    write(11,*) pul_dl/2d0, pul_dl/2d0
    write(11,*)-pul_dl/2d0, pul_dl/2d0
    write(11,*)-pul_dl/2d0,-pul_dl/2d0
    
    if (pul_include_TLM_return) then
    
      write(11,*)
      write(11,*)
      
      wire=pul_nwires
    
! wire center coordinates    
      x=pul_wire_spec(wire)%xc
      y=pul_wire_spec(wire)%yc
      r=pul_wire_spec(wire)%rw

! loop over theta    

      do tloop=0,50
    
        t=2d0*pi*dble(tloop)/50d0
        xp=x+r*cos(t)
        yp=y+r*sin(t)
        write(11,*)xp,yp

      end do

    end if
    
    close(unit=11)  
  
    open(unit=11,file='pul_LC_wires.dat')  
    
    do wire=1,pul_nwires_in
    
! wire center coordinates    
      x=pul_wire_spec(wire)%xc
      y=pul_wire_spec(wire)%yc
      r=pul_wire_spec(wire)%rw

! loop over theta    

      do tloop=0,50
    
        t=2d0*pi*dble(tloop)/50d0
        xp=x+r*cos(t)
        yp=y+r*sin(t)
        write(11,*)xp,yp

      end do
      write(11,*)
      write(11,*)      
    
    end do ! next wire
     
    close(unit=11)  
  
    open(unit=11,file='pul_LC_insulation.dat')  
    
    do wire=1,pul_nwires_in
    
! wire center coordinates    
      x=pul_wire_spec(wire)%xc
      y=pul_wire_spec(wire)%yc
      r=pul_wire_spec(wire)%ri

! loop over theta    

      do tloop=0,50
    
        t=2d0*pi*dble(tloop)/50d0
        xp=x+r*cos(t)
        yp=y+r*sin(t)
        write(11,*)xp,yp

      end do
      write(11,*)
      write(11,*)      
    
    end do ! next wire
     
    close(unit=11)  
    
  end if
  
  return
  
  end SUBROUTINE write_geometry
  
