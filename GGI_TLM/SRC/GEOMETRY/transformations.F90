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
! Name
!     apply_transformation
!
! Description
!     apply transformation to coordinates
!
! Comments:
!      
!
! History
!
!     started 12/02/09 CJS
!

SUBROUTINE apply_transformation(point,trans)

USE geometry_types

! variables passed to subroutine
type(xyz) :: point
type(transformation_type) :: trans

! local variables

real*8 tx,ty,tz
type(xyz) :: translation_point

real*8 		:: mat(3,3)
type(xyz) :: point2


! START

  tx=trans%parameters(1)
  ty=trans%parameters(2)
  tz=trans%parameters(3)

  translation_point%x=trans%parameters(4)
  translation_point%y=trans%parameters(5)
  translation_point%z=trans%parameters(6)

! get rotation matrix
  call euler_angles_to_rotation_matrix(tx,ty,tz,mat)

! rotate first
  call rot(point,mat,point2)

! now translate
  call translate(point2,translation_point,point)

  RETURN

END SUBROUTINE apply_transformation

!
! Name
!     apply_rotation_only
!
! Description
!     apply rotation_only to coordinates
!
! Comments:
!      
!
! History
!
!     started 7/04/09 CJS
!

SUBROUTINE apply_rotation_only(point,trans)

USE geometry_types

! variables passed to subroutine
type(xyz) :: point
type(transformation_type) :: trans

! local variables

real*8 tx,ty,tz

real*8 mat(3,3)
type(xyz) :: point2

! START

  tx=trans%parameters(1)
  ty=trans%parameters(2)
  tz=trans%parameters(3)

! get rotation matrix
  call euler_angles_to_rotation_matrix(tx,ty,tz,mat)

! rotate first
  call rot(point,mat,point2)

  point%x=point2%x
  point%y=point2%y
  point%z=point2%z

RETURN

END SUBROUTINE apply_rotation_only

!
! Name
!     apply_reverse_transformation
!
! Description
!     apply transformation to coordinates
!
! Comments:
!      
!
! History
!
!     started 13/02/09 CJS
!

subroutine apply_reverse_transformation(point,trans)

USE geometry_types

IMPLICIT NONE

! variables passed to subroutine
type(xyz) :: point
type(transformation_type) :: trans

! local variables

real*8 tx,ty,tz
type(xyz) :: translation_point

real*8 mat(3,3)
type(xyz) :: point2

! START

  tx=trans%parameters(1)
  ty=trans%parameters(2)
  tz=trans%parameters(3)

  translation_point%x=-trans%parameters(4)
  translation_point%y=-trans%parameters(5)
  translation_point%z=-trans%parameters(6)

! get rotation matrix
  call euler_angles_to_rotation_matrix(tx,ty,tz,mat)
  
! reverse translate
  call translate(point2,translation_point,point)

! reverse rotate 
  call rotr(point2,mat,point)

RETURN

END SUBROUTINE apply_reverse_transformation
!
!
!
SUBROUTINE rot(point,mat,point2)

USE geometry_types

IMPLICIT NONE

! rotate coordinate system

type(xyz) :: point
real*8 mat(3,3)
type(xyz) :: point2

point2%x=point%x*mat(1,1)+point%y*mat(1,2)+point%z*mat(1,3)
point2%y=point%x*mat(2,1)+point%y*mat(2,2)+point%z*mat(2,3)
point2%z=point%x*mat(3,1)+point%y*mat(3,2)+point%z*mat(3,3)

RETURN

END SUBROUTINE rot
!
!
!
SUBROUTINE rotr(point,mat,point2)

USE geometry_types

IMPLICIT NONE

! reverse rotation of coordinate system

type(xyz) :: point
real*8 mat(3,3)
type(xyz) :: point2

! START

point2%x=point%x*mat(1,1)+point%y*mat(2,1)+point%z*mat(3,1)
point2%y=point%x*mat(1,2)+point%y*mat(2,2)+point%z*mat(3,2)
point2%z=point%x*mat(1,3)+point%y*mat(2,3)+point%z*mat(3,3)

RETURN

END SUBROUTINE rotr
!
!
!
SUBROUTINE euler_angles_to_rotation_matrix(tx,ty,tz,mat)

USE geometry_types

IMPLICIT NONE

! calculate rotation matrix from euler angles

real*8 tx,ty,tz
real*8 mat(3,3)

real*8 cx,cy,cz,sx,sy,sz

! START

cx=cos(tx)
cy=cos(ty)
cz=cos(tz)

sx=sin(tx)
sy=sin(ty)
sz=sin(tz)

mat(1,1)=cy*cz
mat(1,2)=cy*sz
mat(1,3)=-sy

mat(2,1)=sx*sy*cz-cx*sz
mat(2,2)=sx*sy*sz+cx*cz
mat(2,3)=cy*sx

mat(3,1)=cx*sy*cz+sx*sz
mat(3,2)=cx*sy*sz-sx*cz
mat(3,3)=cy*cx

RETURN

END SUBROUTINE euler_angles_to_rotation_matrix
!
!
!
SUBROUTINE translate(point,translation_point,point2)

USE geometry_types

IMPLICIT NONE

! translate coordinate
type(xyz) :: point
type(xyz) :: translation_point
type(xyz) :: point2

!START

point2%x=point%x+translation_point%x
point2%y=point%y+translation_point%y
point2%z=point%z+translation_point%z

RETURN

END SUBROUTINE translate
!
!
!
SUBROUTINE xyz_point_to_rthetaphi(point,r,theta,phi)

USE geometry_types
USE constants

IMPLICIT NONE

! calculate rotation matrix from euler angles

type(xyz)	:: point
real*8 		:: r,theta,phi

! START

  r=sqrt(point%x*point%x+point%y*point%y+point%z*point%z)
  if (point%x.ne.0d0) then
    phi=atan2(point%y,point%x)
  else
    if (point%y.gt.0d0) then
      phi=pi/2d0
    else if (point%y.lt.0d0) then
      phi=-pi/2d0
    else
      phi=0d0
    end if
  end if
  if(r.ne.0d0) then
    theta=acos(point%z/r)
  else
    theta=0d0
  end if
  
  return
  
end
!
!
!
subroutine rthetaphi_to_xyz_point(r,theta,phi,point)

USE geometry_types
USE constants

IMPLICIT NONE

real*8 r,theta,phi

type(xyz)	:: point

! START

  point%x=r*sin(theta)*cos(phi)
  point%y=r*sin(theta)*sin(phi)
  point%z=r*cos(theta)
  
  return
  
end    
!
!
!
subroutine rphiz_to_xyz_point(r,phi,z,point)

USE geometry_types
USE constants

IMPLICIT NONE

real*8 r,phi,z

type(xyz)	:: point

! START

  point%x=r*cos(phi)
  point%y=r*sin(phi)
  point%z=z
  
  return
  
end    
