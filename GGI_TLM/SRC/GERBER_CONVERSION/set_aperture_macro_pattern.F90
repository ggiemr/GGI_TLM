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
!     SUBROUTINE set_aperture_macro_pattern
!     SUBROUTINE reverse_transform
!
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 8/2/19 CJS
!     
!
SUBROUTINE set_aperture_macro_pattern(rloop)

USE gerber
USE evaluate_string_expression
!
IMPLICIT NONE

integer :: rloop

! local variables

integer iax,iay
real*8  :: lx,ly,rx,ry
real*8  :: ax,ay,ar
integer :: ns
real*8  :: theta,theta_0,alpha

real*8 :: xp,yp,xn,yn,pline,ptest
real*8  :: value
integer :: level
logical :: verbose_expression=.FALSE.   ! used in evaluating expressions
integer :: i

integer :: AM    ! aperture macro number
integer :: np,primitive  ! number of aperture macro primitives, primitive loop variable
integer :: lp            ! position in list of primitives
integer :: ptype         ! primitive type
integer :: nm            ! number of modifiers

real*8  :: modifier_values(1:maxAM_modifiers)
real*8  :: primitive_size

integer :: exposure
real*8  :: radius,cx,cy,rot,x1,y1,x2,y2,w,h
real*8  :: ldx,ldy,llen,lnx,lny,dline

integer :: var_number
real*8  :: var_value

logical :: set_point

integer :: loop

real*8,parameter :: pi=3.1415926535d0

! set the aperture pattern

! START

! Set any variables defined for this aperture

write(*,*)'Setting aperture from aperture macro'
write(*,*)'Aperture=',aperture,' n_vars=',A_nvars(aperture)
write(*,*)' var     value '

n_vars=A_nvars(aperture)
ALLOCATE( var_list(1:n_vars) ) ! note var_list is for numerical values

do loop=1,rloop   ! the fisrt loop works out the size of the aperture and the second allocates and sets the pattern

  write(*,*)
  if (loop.EQ.1) then
    write(*,*)'Loop=1, working out the approximate aperture size'
  else
    write(*,*)'Loop=2, setting aperture'  
  end if

  do i=1,n_vars
    level=0
    CALL eval_expression(Avars(aperture,i),value,level,verbose)
    var_list(i)=value
    write(*,*)i,var_list(i)
  end do

  AM=A_AMnumber(aperture)
  np=AM_n_primitives(AM)

  write(*,*)'Aperture Macro number is:',AM
  write(*,*)'Aperture Macro name is  :',trim(AMname(AM))
  
  if (loop.EQ.1) then 
! work out aperture size starting from an initial value of zero
    Asize(aperture)=0d0
    
  else
! we have the aperture size so allocate memory for the pattern

    if (allocated( ap )) Deallocate( ap )

    anx=NINT( 0.5d0*Asize(aperture)/dl )+1
    any=NINT( 0.5d0*Asize(aperture)/dl )+1
 
    write(*,*)'Aperture size is:',Asize(aperture)
    write(*,*)'dl=',dl
    write(*,*)'anx=',anx,' any=',any
 
    ALLOCATE( ap(-anx:anx,-any:any) )
    ap(-anx:anx,-any:any)=0
 
  end if ! loop.EQ.1
  
! loop over the primitives which define the aperture
  
  do primitive=1,np
  
    lp=AM_to_primitive_list(AM,primitive)
    ptype=AM_primitive_number(lp)
    nm=n_AM_modifiers(lp)
    
    write(*,*)'Aperture macro primitive number:',primitive
    write(*,*)'position in primitive list     :',lp
    write(*,*)'primitive type                 :',ptype
    write(*,*)'number of modifiers            :',nm
    
    if (ptype.LT.0) then
! this 'primitive' is actually a variable definition

      var_number=-ptype
      write(*,*)'Reset variable value for variable number',var_number
     
      CALL eval_expression(AM_modifiers(lp,1),var_value,level,verbose_expression)
      write(*,*)'Variable value=',var_value
     
      if (var_number.GT.n_vars) then
        write(*,*)'ERROR: var_number is greater than n_vars, n_vars=',n_vars
        STOP
      end if
      
      var_list(var_number)=var_value
    
    else if (ptype.EQ.1) then
    
      write(*,*)'Set circle primitive'
! Get the parameters for a circular primitive
      
      if ( (nm.NE.4).AND.(nm.NE.5) ) then
        write(*,*)'ERROR: there should be 4 or 5 modifiers for a circular primitive, found ',nm
        STOP
      end if
! evaluate the modifiers
      write(*,*)'Modifier list:'
      do i=1,nm
      
        level=0
        CALL eval_expression(AM_modifiers(lp,i),modifier_values(i),level,verbose_expression)
        write(*,*)i,trim(AM_modifiers(lp,i)),modifier_values(i)
      
      end do 
      
      exposure=NINT(modifier_values(1))
      radius=(modifier_values(2)/2d0)*s
      cx=modifier_values(3)*s
      cy=modifier_values(4)*s
      if (nm.EQ.4) then
        rot=0d0
      else
        rot=modifier_values(5)
      end if
      
      if (loop.EQ.1) then
      
        primitive_size=2d0*(sqrt(cx**2+cy**2)+radius)
        Asize(aperture)=max(primitive_size,Asize(aperture))

      else
      
! loop over the aperture points setting the pattern    
        do iax=-anx,anx
          do iay=-any,any  
! get the coordinates of this point on the aperture
  
            ax=iax*dl
            ay=iay*dl
          
            CALL reverse_transform(ax,ay,cx,cy,rot)  ! reverse the translation and rotation
    
            ar=sqrt(ax*ax+ay*ay)
        
            if  (ar.LE.radius) then
              ap(iax,iay)=exposure
            end if
        
          end do
        end do
      end if    ! loop
      
    else if (ptype.EQ.20) then
      write(*,*)'Set vector line primitive'
      
! Get the parameters for a vector line primitive
      
      if (nm.NE.7) then
        write(*,*)'ERROR: there should be 7 modifiers for a vector line primitive, found ',nm
        STOP
      end if
! evaluate the modifiers
      write(*,*)'Modifier list:'
      do i=1,nm
      
        level=0
        CALL eval_expression(AM_modifiers(lp,i),modifier_values(i),level,verbose_expression)
        write(*,*)i,trim(AM_modifiers(lp,i)),modifier_values(i)
      
      end do 
      
      exposure=NINT(modifier_values(1))
      w=modifier_values(2)*s
      x1=modifier_values(3)*s
      y1=modifier_values(4)*s
      x2=modifier_values(5)*s
      y2=modifier_values(6)*s
      rot=modifier_values(7)
      
      ldx=x2-x1
      ldy=y2-y1
      llen=sqrt(ldx*ldx+ldy*ldy)
      ldx=ldx/llen
      ldy=ldy/llen
      lnx=-ldy
      lny=ldx
      cx=0d0
      cy=0d0
      
      if (loop.EQ.1) then
      
        primitive_size=2d0*(max(sqrt(x1**2+y1**2),sqrt(x2**2+y2**2))+w)
        Asize(aperture)=max(primitive_size,Asize(aperture))

      else

! loop over the aperture points  setting the pattern      
        do iax=-anx,anx
          do iay=-any,any  
! get the coordinates of this point on the aperture
  
            ax=iax*dl
            ay=iay*dl
          
            CALL reverse_transform(ax,ay,cx,cy,rot)  ! reverse the translation and rotation

! calculate distance from point to a line from (x1,y1) to (x2,y2)           
            dline=ax*lnx+ay*lny
! calculate projection of vector from (x1,y1) to point, along the line
            h=(ax-x1)*ldx+(ay-y1)*ldy
            
            if ( (abs(dline).LE.w/2d0).AND.       &
                 (h.GE.0d0).AND.(h.LT.llen) ) then
              ap(iax,iay)=exposure
            end if
        
          end do
        end do
      end if    ! loop
    
    else if (ptype.EQ.21) then
      write(*,*)'Set centre line primitive'
! Get the parameters for a centre line primitive
      
      if (nm.NE.6) then
        write(*,*)'ERROR: there should be 6 modifiers for a centre line primitive, found ',nm
        STOP
      end if
! evaluate the modifiers
      write(*,*)'Modifier list:'
      do i=1,nm
      
        level=0
        CALL eval_expression(AM_modifiers(lp,i),modifier_values(i),level,verbose_expression)
        write(*,*)i,trim(AM_modifiers(lp,i)),modifier_values(i)
      
      end do 
      
      exposure=NINT(modifier_values(1))
      w=modifier_values(2)*s
      h=modifier_values(3)*s
      cx=modifier_values(4)*s
      cy=modifier_values(5)*s
      rot=modifier_values(6)
      
      if (loop.EQ.1) then
      
        primitive_size=2d0*(sqrt(cx**2+cy**2)+sqrt(w*w+h*h)/2d0)
        Asize(aperture)=max(primitive_size,Asize(aperture))

      else

! loop over the aperture points  setting the pattern      
        do iax=-anx,anx
          do iay=-any,any  
! get the coordinates of this point on the aperture
  
            ax=iax*dl
            ay=iay*dl
          
            CALL reverse_transform(ax,ay,cx,cy,rot)  ! reverse the translation and rotation
           
            if ( (abs(ax).LE.w/2d0).AND. &
                 (abs(ay).LE.h/2d0)       ) then
              ap(iax,iay)=exposure
            end if
        
          end do
        end do
      end if    ! loop
    
    else if (ptype.EQ.4) then
      write(*,*)'Set outline primitive'  
          
        AMflag=.TRUE.
        AMflag_outline=.TRUE.

    else if (ptype.EQ.5) then
      write(*,*)'Set polygon primitive'   
! Get the parameters for a polygon primitive
      
      if (nm.NE.6) then
        write(*,*)'ERROR: there should be 6 modifiers for a polygon primitive, found ',nm
        STOP
      end if
! evaluate the modifiers
      write(*,*)'Modifier list:'
      do i=1,nm
      
        level=0
        CALL eval_expression(AM_modifiers(lp,i),modifier_values(i),level,verbose_expression)
        write(*,*)i,trim(AM_modifiers(lp,i)),modifier_values(i)
      
      end do 
      
      exposure=NINT(modifier_values(1))
      ns=NINT(modifier_values(2))
      cx=modifier_values(3)*s
      cy=modifier_values(4)*s
      radius=(modifier_values(5)/2d0)*s
      rot=modifier_values(6)
      
      if (loop.EQ.1) then
      
        primitive_size=2d0*(sqrt(cx**2+cy**2)+radius)
        Asize(aperture)=max(primitive_size,Asize(aperture))

      else
      
! loop over the aperture points setting the pattern    
        do iax=-anx,anx
          do iay=-any,any  
! get the coordinates of this point on the aperture
  
            ax=iax*dl
            ay=iay*dl
          
            CALL reverse_transform(ax,ay,cx,cy,rot)  ! reverse the translation and rotation
    
            set_point=.FALSE.
            ar=sqrt(ax*ax+ay*ay)
        
            if  (ar.LE.radius) then
              set_point=.TRUE.
            end if
            
! loop over sides
            
            do i=1,ns
        
              theta=(dble(i)-0.5d0)*2d0*pi/dble(ns) ! angle to the midpoint of the edge
              xn=cos(theta)  ! normal to the edge
              yn=sin(theta)
          
              alpha=dble(i-1)*2d0*pi/dble(ns)  ! angle to a point on the edge (first vertex)
          
              xp=radius*cos(alpha)  ! get a point on the edge (first vertex)
              yp=radius*sin(alpha)
          
              pline=xp*xn+yp*yn   ! point on edge dotted with normal direction
          
              ptest=ax*xn+ay*yn   ! point dotted with normal direction
              if (ptest.GT.pline) set_point=.FALSE.
        
            end do ! next side
        
            if (set_point) then
              ap(iax,iay)=exposure
            end if
        
          end do ! next ix      
        end do ! next iy
        
      end if    ! loop
          
    else if (ptype.EQ.6) then
      write(*,*)'Set Moire primitive'
          
        AMflag=.TRUE.
        AMflag_Moire=.TRUE.
    
    else if (ptype.EQ.7) then
      write(*,*)'Set thermal primitive'
          
        AMflag=.TRUE.
        AMflag_thermal=.TRUE.
   
    else if (ptype.NE.0) then
    
      write(*,*)'ERROR: unknown primitive type:',ptype
      STOP
      
    end if    
    
  end do ! next primitive in this aperture macro
  
end do ! loop

! finish up

DEALLOCATE( var_list )

RETURN
          
END SUBROUTINE set_aperture_macro_pattern
!
! ____________________________________________________________________________
!
!
SUBROUTINE reverse_transform(x,y,cx,cy,rot)  ! reverse the translation and rotation

real*8 :: x,y,cx,cy,rot

! local variables

real*8 :: xt,yt
real*8 :: rad_rot
real*8,parameter :: pi=3.1415926535d0

! START

! reverse the rotation

theta=-rot*pi/180d0     ! reverse the angle and convert from degrees to radians

xt=cos(theta)*x-sin(theta)*y
yt=sin(theta)*x+cos(theta)*y

! remove the offset

x=xt-cx
y=yt-cy

RETURN

end SUBROUTINE reverse_transform
