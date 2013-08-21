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
!SUBROUTINE is_point_inside_tet
!SUBROUTINE line_triangle_intersection
!SUBROUTINE triangle_normal
!SUBROUTINE triangle_unit_normal
!SUBROUTINE vector_product
!
!FUNCTION   dot
!FUNCTION   vector_length
!FUNCTION   distance
!FUNCTION   xyz_dot
!FUNCTION   xyz_vector_length
!FUNCTION   xyz_distance
!
!
SUBROUTINE is_point_inside_tet(point,tet,inside)
!
!
!
USE geometry_types
USE geometry_operators
USE constants

IMPLICIT NONE

  type(xyz)		:: point
  type(xyz_tet)		:: tet
  logical		:: inside

! local variables

  type(xyz)		:: v1,v2,v3,vp
  type(xyz)		:: normal
  real*8		:: test1,test2

! function variables

  real*8		:: xyz_dot
  real*8		:: xyz_vector_length
  type(xyz)		:: xyz_vector_product

! START

!  write(*,*)'CALLED:is_point_inside_tet'
!  write(*,*)'point',point%x,point%y,point%z
!  write(*,*)'V1',tet%vertex(1)%x,tet%vertex(1)%y,tet%vertex(1)%z
!  write(*,*)'V2',tet%vertex(2)%x,tet%vertex(2)%y,tet%vertex(2)%z
!  write(*,*)'V3',tet%vertex(3)%x,tet%vertex(3)%y,tet%vertex(3)%z
!  write(*,*)'V4',tet%vertex(4)%x,tet%vertex(4)%y,tet%vertex(4)%z

! assume point is outside
  inside=.FALSE.
  
! test 1 - are tet point4 and the test point the same side of the plane containing points 1,2 and 3?
  v1=tet%vertex(2)-tet%vertex(1)
  v2=tet%vertex(3)-tet%vertex(1)
  normal=xyz_vector_product(v1,v2)
  v3=tet%vertex(4)-tet%vertex(1)
  vp=point-tet%vertex(1)
  test1=xyz_dot(normal,v3)
  test2=xyz_dot(normal,vp)
  
!  write(*,*)'test 1',test1*test2

! if the signs are different then the test fails so return fail
  if (test1*test2.lt.0d0) RETURN
  
! test 2 - are tet point3 and the test point the same side of the plane containing points 1,2 and 4?
  v1=tet%vertex(2)-tet%vertex(1)
  v2=tet%vertex(4)-tet%vertex(1)
  normal=xyz_vector_product(v1,v2)
  v3=tet%vertex(3)-tet%vertex(1)
  vp=point-tet%vertex(1)
  test1=xyz_dot(normal,v3)
  test2=xyz_dot(normal,vp)
  
!  write(*,*)'test 2',test1*test2

! if the signs are different then the test fails so return fail
  if (test1*test2.lt.0d0) RETURN
  
! test 3 - are tet point2 and the test point the same side of the plane containing points 1,3 and 4?
  v1=tet%vertex(3)-tet%vertex(1)
  v2=tet%vertex(4)-tet%vertex(1)
  normal=xyz_vector_product(v1,v2)
  v3=tet%vertex(2)-tet%vertex(1)
  vp=point-tet%vertex(1)
  test1=xyz_dot(normal,v3)
  test2=xyz_dot(normal,vp)
  
!  write(*,*)'test 3',test1*test2

! if the signs are different then the test fails so return fail
  if (test1*test2.lt.0d0) RETURN
  
! test 4 - are tet point1 and the test point the same side of the plane containing points 2,3 and 4?
  v1=tet%vertex(3)-tet%vertex(2)
  v2=tet%vertex(4)-tet%vertex(2)
  normal=xyz_vector_product(v1,v2)
  v3=tet%vertex(1)-tet%vertex(2)
  vp=point-tet%vertex(2)
  test1=xyz_dot(normal,v3)
  test2=xyz_dot(normal,vp)
  
!  write(*,*)'test 4',test1*test2

! if the signs are different then the test fails so return fail
  if (test1*test2.lt.0d0) RETURN
  
  inside=.TRUE.

  RETURN

END SUBROUTINE is_point_inside_tet
!
!
!
SUBROUTINE line_triangle_intersection(line_point1,line_point2,triangle,	&
			              intersection_point,intersection)
!
! The algorithm is described on the web page:
!				      
!http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle%28%29
!
USE geometry_types
USE geometry_operators
USE constants

IMPLICIT NONE

  type(xyz)		:: line_point1,line_point2
  type(xyz_triangle)	:: triangle
  type(xyz)		:: intersection_point
  logical		:: intersection

! local variables

  type(xyz)		:: normal

  type(xyz)		:: V1_P1
  type(xyz)		:: P2_P1
  
  type(xyz)		:: u,v,w
  
  real*8		:: r
  real*8		:: normal_length
  real*8		:: num,den
  
  real*8		:: uv,wv,vv,wu,uu
  real*8		:: d
  real*8		:: s,t

! function variables

  real*8		:: xyz_dot
  real*8		:: xyz_vector_length
  type(xyz)		:: xyz_vector_product

! START

! assume no intersection
  intersection=.FALSE.

! get the normal to the triangle
  CALL triangle_normal(triangle,normal)
  
  normal_length=xyz_vector_length(normal)
  
! Triangle has zero area so return false
  if (normal_length.lt.small) RETURN

! calculate unit normal  
  normal%x=normal%x/normal_length
  normal%y=normal%y/normal_length
  normal%z=normal%z/normal_length

! check that the line from point1 to point2 intersects the plane of the triangle

  V1_P1=triangle%vertex(1)-line_point1
  P2_P1=line_point2       -line_point1
  
  num=xyz_dot(normal,V1_P1)
  den=xyz_dot(normal,P2_P1)
  
! line is parallel to the plane so return false
  if (den.eq.0d0) RETURN
  
  r=num/den
  
  intersection_point%x=line_point1%x+r*(line_point2%x-line_point1%x)
  intersection_point%y=line_point1%y+r*(line_point2%y-line_point1%y)
  intersection_point%z=line_point1%z+r*(line_point2%z-line_point1%z)
  
  u=triangle%vertex(2)-triangle%vertex(1)
  v=triangle%vertex(3)-triangle%vertex(1)
  w=intersection_point-triangle%vertex(1)
  
  uv=xyz_dot(u,v)
  wv=xyz_dot(w,v)
  vv=xyz_dot(v,v)
  wu=xyz_dot(w,u)
  uu=xyz_dot(u,u)
  
  d=(uv*uv-uu*vv)
  
  s=(uv*wv-vv*wu)/d
  
!  if ((s.LT.0d0).OR.(s.GT.1d0) ) RETURN  ! exact check - changed to include some tolerance
  if ((s.LT.-small).OR.(s.GT.1d0+small) ) RETURN
  
  t=(uv*wu-uu*wv)/d
  
!  if ( (t.LT.0d0).OR.(s+t.GT.1d0) ) RETURN  ! exact check - changed to include some tolerance
  if ( (t.LT.-small).OR.(s+t.GT.1d0+small) ) RETURN
   
  intersection=.TRUE.
 
  RETURN

END SUBROUTINE line_triangle_intersection
!
!
!

SUBROUTINE triangle_normal(triangle,normal)

USE geometry_types

IMPLICIT NONE

  type(xyz_triangle)	:: triangle
  type(xyz)		:: normal
  
! local variables

  real*8		:: vx1,vy1,vz1
  real*8		:: vx2,vy2,vz2
  real*8		:: xn,yn,zn

! function variables

  real*8		:: vector_length

! START
  
  vx1=triangle%vertex(3)%x-triangle%vertex(1)%x
  vy1=triangle%vertex(3)%y-triangle%vertex(1)%y
  vz1=triangle%vertex(3)%z-triangle%vertex(1)%z

  vx2=triangle%vertex(2)%x-triangle%vertex(1)%x
  vy2=triangle%vertex(2)%y-triangle%vertex(1)%y
  vz2=triangle%vertex(2)%z-triangle%vertex(1)%z
  
  CALL vector_product(vx1,vy1,vz1,vx2,vy2,vz2,xn,yn,zn)
  
  normal%x=xn
  normal%y=yn
  normal%z=zn
  
  RETURN
  
END SUBROUTINE triangle_normal
!
!
!

SUBROUTINE triangle_unit_normal(triangle,normal)

USE geometry_types

IMPLICIT NONE

  type(xyz_triangle)	:: triangle
  type(xyz)		:: normal
  
! local variables

  real*8		:: vx1,vy1,vz1
  real*8		:: vx2,vy2,vz2
  real*8		:: xn,yn,zn
  
  real*8		:: norm

! function variables

  real*8		:: vector_length

! START

  
  vx1=triangle%vertex(3)%x-triangle%vertex(1)%x
  vy1=triangle%vertex(3)%y-triangle%vertex(1)%y
  vz1=triangle%vertex(3)%z-triangle%vertex(1)%z

  vx2=triangle%vertex(2)%x-triangle%vertex(1)%x
  vy2=triangle%vertex(2)%y-triangle%vertex(1)%y
  vz2=triangle%vertex(2)%z-triangle%vertex(1)%z
  
  CALL vector_product(vx1,vy1,vz1,vx2,vy2,vz2,xn,yn,zn)
  
  norm=vector_length(xn,yn,zn)
  
  if (norm.ne.0d0) then
  
    normal%x=xn/norm
    normal%y=yn/norm
    normal%z=zn/norm
    
  else
  
    GOTO 9000
  
  end if
  
  RETURN
  
9000  CALL write_line('ERROR in triangle_normal',0,.TRUE.)
      CALL write_line('triangle area =0',0,.TRUE.)
      STOP

  
END SUBROUTINE triangle_unit_normal
!
!
!

SUBROUTINE vector_product(ax,ay,az,bx,by,bz,cx,cy,cz)


USE constants

IMPLICIT NONE

real*8 ax,ay,az,bx,by,bz,cx,cy,cz

! START

  cx=ay*bz-az*by
  cy=az*bx-ax*bz
  cz=ax*by-ay*bx
  
  RETURN
  
END    
!
!
!
 FUNCTION dot(ax,ay,az,bx,by,bz) RESULT(res)

IMPLICIT NONE

 real*8 ax,ay,az,bx,by,bz
 real*8 res
 
! START
 
 res=ax*bx+ay*by+az*bz
 
 RETURN
 END
!
!
!
 FUNCTION vector_length(ax,ay,az) RESULT(res)

IMPLICIT NONE

 real*8 ax,ay,az
 real*8 res
 
! START
 
 res=sqrt(ax*ax+ay*ay+az*az)
 
 return
 end

FUNCTION xyz_vector_product(a,b) RESULT(res)

USE geometry_types

IMPLICIT NONE
  type(xyz)		:: a,b
  type(xyz)		:: res

! START

  res%x=a%y*b%z-a%z*b%y
  res%y=a%z*b%x-a%x*b%z
  res%z=a%x*b%y-a%y*b%x
  
  RETURN
  
END    
!
!
!
FUNCTION xyz_dot(a,b) RESULT(res)

USE geometry_types

IMPLICIT NONE

  type(xyz)		:: a,b
  real*8 		:: res
 
! START
 
 res=a%x*b%x+a%y*b%y+a%z*b%z
 
 
 RETURN
 END
!
!
!
 FUNCTION xyz_vector_length(a) RESULT(res)

USE geometry_types

IMPLICIT NONE

  type(xyz)		:: a
  real*8		:: res
 
! START
 
 res=sqrt(a%x*a%x+a%y*a%y+a%z*a%z)
 
 RETURN
 END
!
!
!
 FUNCTION xyz_distance(a,b) RESULT(res)

USE geometry_types
USE geometry_operators

IMPLICIT NONE

  type(xyz)		:: a,b
  real*8		:: res

! local variables  
  type(xyz)		:: c
  
! function variables  
  real*8		:: xyz_vector_length
 
! START

 c=a-b
 
 res=xyz_vector_length(c)
 
 RETURN
 END
