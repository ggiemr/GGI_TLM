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
!SUBROUTINE write_tube_points_vtk
!SUBROUTINE write_cone_points_vtk
!SUBROUTINE get_LT_vectors
!
! NAME
!     SUBROUTINE write_tube_points_vtk
!
! DESCRIPTION
!     write_tube_points_vtk
!
!     write_tube_points vtk format file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 31/08/2012 CJS
!
!
SUBROUTINE write_tube_points_vtk(x1,y1,z1,x2,y2,z2,offset,number_of_surfaces,file_unit)

USE constants

IMPLICIT NONE

real*8		:: x1,y1,z1,x2,y2,z2
real*8		:: offset
integer 	:: number_of_surfaces
integer		:: file_unit

! local variables

real*8		:: vx,vy,vz
real*8		:: vn1x,vn1y,vn1z
real*8		:: vn2x,vn2y,vn2z

integer		:: n,n_theta
real*8		:: theta_min,theta_max,d_theta,theta1,theta2

! START
  
! get vector along line and two orthogonal transverse vectors
    CALL get_LT_vectors(x1,y1,z1,x2,y2,z2,vx,vy,vz,vn1x,vn1y,vn1z,vn2x,vn2y,vn2z)
  
    theta_min=0d0
    theta_max=2d0*pi
  
    n_theta=number_of_surfaces

! loop over theta

    d_theta=(theta_max-theta_min)/dble(n_theta)

    do n=1,n_theta  
  
      theta1=dble(n-1)*d_theta
      theta2=dble(n)*d_theta    
    
      write(file_unit,8000)	x1+offset*(vn1x*cos(theta1)+vn2x*sin(theta1)),	&
    				y1+offset*(vn1y*cos(theta1)+vn2y*sin(theta1)),	&
    				z1+offset*(vn1z*cos(theta1)+vn2z*sin(theta1))
    
      write(file_unit,8000)	x1+offset*(vn1x*cos(theta2)+vn2x*sin(theta2)),	&
    				y1+offset*(vn1y*cos(theta2)+vn2y*sin(theta2)),	&
    				z1+offset*(vn1z*cos(theta2)+vn2z*sin(theta2))
    
      write(file_unit,8000)	x2+offset*(vn1x*cos(theta2)+vn2x*sin(theta2)),   &
    				y2+offset*(vn1y*cos(theta2)+vn2y*sin(theta2)),	&
    				z2+offset*(vn1z*cos(theta2)+vn2z*sin(theta2))
    
      write(file_unit,8000)	x2+offset*(vn1x*cos(theta1)+vn2x*sin(theta1)),	&
    				y2+offset*(vn1y*cos(theta1)+vn2y*sin(theta1)),	&
    				z2+offset*(vn1z*cos(theta1)+vn2z*sin(theta1))
    
    end do
	      
8000  format(3E14.5)
  
  RETURN
  
  
END SUBROUTINE write_tube_points_vtk
!
! NAME
!     SUBROUTINE write_cone_points_vtk
!
! DESCRIPTION
!     write_cone_points_vtk
!
!     write_cone_points vtk format file
!     Used to plot vector fields
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/02/2013 CJS
!
!
SUBROUTINE write_cone_points_vtk(x1,y1,z1,x2,y2,z2,offset,number_of_surfaces,file_unit)

USE constants

IMPLICIT NONE

real*8		:: x1,y1,z1,x2,y2,z2
real*8		:: offset
integer 	:: number_of_surfaces
integer		:: file_unit

! local variables

real*8		:: vx,vy,vz
real*8		:: vn1x,vn1y,vn1z
real*8		:: vn2x,vn2y,vn2z

integer		:: n,n_theta
real*8		:: theta_min,theta_max,d_theta,theta1,theta2

! START
  
! get vector along line and two orthogonal transverse vectors
    CALL get_LT_vectors(x1,y1,z1,x2,y2,z2,vx,vy,vz,vn1x,vn1y,vn1z,vn2x,vn2y,vn2z)
  
    theta_min=0d0
    theta_max=2d0*pi
  
    n_theta=number_of_surfaces

! loop over theta

    d_theta=(theta_max-theta_min)/dble(n_theta)

    do n=1,n_theta  
  
      theta1=dble(n-1)*d_theta
      theta2=dble(n)*d_theta    
    
     write(file_unit,8000)	x1+offset*(vn1x*cos(theta2)+vn2x*sin(theta2)),	&
    				y1+offset*(vn1y*cos(theta2)+vn2y*sin(theta2)),	&
    				z1+offset*(vn1z*cos(theta2)+vn2z*sin(theta2))
    
      write(file_unit,8000)	x2,  &
    				y2,  &
    				z2
    
      write(file_unit,8000)	x2,  &
    				y2,  &
    				z2
    
      write(file_unit,8000)	x1+offset*(vn1x*cos(theta1)+vn2x*sin(theta1)),	&
    				y1+offset*(vn1y*cos(theta1)+vn2y*sin(theta1)),	&
    				z1+offset*(vn1z*cos(theta1)+vn2z*sin(theta1))
    
     end do
	      
8000  format(3E14.5)
  
  RETURN
  
  
END SUBROUTINE write_cone_points_vtk
!
!
! Name 
!     
!
! Description
!     return unit vectors along the line between two points and two mutually
!     orthogonal transverse vectors
!
! Comments:
!      
!
! History
!
!     started 11/08/10 CJS
!

  SUBROUTINE get_LT_vectors(xmin,ymin,zmin,xmax,ymax,zmax,vx,vy,vz,	&
                      ox1,oy1,oz1,ox2,oy2,oz2)
IMPLICIT NONE

! variables passed to subroutine

  real*8 xmin,ymin,zmin,xmax,ymax,zmax
  real*8 vx,vy,vz,ox1,oy1,oz1,ox2,oy2,oz2

! local_variables
  
  real*8 tx,ty,tz,l1,l2,norm

! START

! vector along wire
  vx=xmax-xmin
  vy=ymax-ymin
  vz=zmax-zmin
  
  norm=sqrt(vx*vx+vy*vy+vz*vz)
  if (norm.eq.0d0) then
    write(*,*)'Error in get LT vectors: Zero length vector'
    stop
  end if
  
  vx=vx/norm
  vy=vy/norm
  vz=vz/norm
  
! calculate normal vector 1.
  tx=1d0
  ty=0d0
  tz=0d0
  
  ox1=vy*tz-vz*ty
  oy1=vz*tx-vx*tz
  oz1=vx*ty-vy*tx
  
  l1=sqrt(ox1*ox1+oy1*oy1+oz1*oz1)
  
! calculate normal vector 2.
  tx=0d0
  ty=1d0
  tz=0d0
  
  ox2=vy*tz-vz*ty
  oy2=vz*tx-vx*tz
  oz2=vx*ty-vy*tx
  
  l2=sqrt(ox2*ox2+oy2*oy2+oz2*oz2)
  
! choose first normal vector as the longest of ox1 and ox2
  if (l2.gt.l1) then
  
    ox1=ox2
    oy1=oy2
    oz1=oz2
    l1=l2
    
  end if
  
  ox1=ox1/l1
  oy1=oy1/l1
  oz1=oz1/l1
      
! calculate second transverse vector as v cross o1  
  tx=ox1
  ty=oy1
  tz=oz1
  
  ox2=vy*tz-vz*ty
  oy2=vz*tx-vx*tz
  oz2=vx*ty-vy*tx
  
  l2=sqrt(ox2*ox2+oy2*oy2+oz2*oz2)
  
  ox2=ox2/l2
  oy2=oy2/l2
  oz2=oz2/l2
  
  return
  
END SUBROUTINE get_LT_vectors
