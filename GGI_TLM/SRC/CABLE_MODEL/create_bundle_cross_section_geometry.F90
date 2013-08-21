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
!SUBROUTINE create_bundle_cross_section_geometry
!SUBROUTINE write_conductor
!SUBROUTINE write_conductor_dash
!SUBROUTINE check_cable_intersection 
!SUBROUTINE calculate_cable_centre_coordinates
!
! NAME
!     SUBROUTINE create_bundle_cross_section
!
! DESCRIPTION
!       Given a list of conductor radii, construct a bundle cross section geometry
!       by adding conductors one at a time whilst trying to minimise the spread of the cross section
!       away from the centre point. This should build a fairly compact cross section.
!   
!       if the bundle radius is greater than max_radius_LC*dl/2 then the return radius is increased.
!       max_radius_LC=0.7 here.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/09/2012 CJS
!
!
SUBROUTINE create_bundle_cross_section_geometry(LC_matrix_dimension,ri,xc,yc)

USE file_information
USE TLM_general

IMPLICIT NONE

integer LC_matrix_dimension,tot_LC_matrix_dimension

real*8 	:: rw(1:LC_matrix_dimension+1)
real*8 	:: ri(1:LC_matrix_dimension+1)
real*8 	:: xc(1:LC_matrix_dimension+1)
real*8 	:: yc(1:LC_matrix_dimension+1)
integer	:: bundle_number

! local variables

integer conductor,n_conductors,tot_n_conductors
integer conductor1,conductor2

real*8 delta

real*8 x0,y0,r0,d0
real*8 x1,y1,r1
real*8 x2,y2,r2
real*8 xmin,ymin,dmin,xmax,ymax,xshift,yshift
logical cable_intersection
logical cable_set

integer i,k
integer side

character ch

real*8 radius,max_radius,min_radius,dielectric_radius_factor,rmax,rmax_required,rmax_found

! START

!  write(*,*)'CALLED: create_bundle_cross_section_geometry',LC_matrix_dimension
  
!  do i=1,LC_matrix_dimension+1
!    write(*,*)i,ri(i),xc(i),yc(i)
!  end do
  
  n_conductors=LC_matrix_dimension
  tot_n_conductors=LC_matrix_dimension+1
  
  delta=0d0
 
! estimate a cable separation 1/20 of the maximum cable radius 
  do i=1,n_conductors
    delta=max(delta,ri(i))
  end do
  delta=delta/20d0
 
  conductor=1

10 CONTINUE

    dmin=1e30
    cable_set=.FALSE.

    r0=ri(conductor)+delta/2d0
    
    if (conductor.eq.1) then
!   put the first conductor at the origin

      xmin=0d0
      ymin=0d0
      cable_set=.TRUE.
  
    else if (conductor.eq.2) then
! put the second conductor touching the first, offset in the x direction

      xmin=ri(conductor-1)+ri(conductor)+delta
      ymin=0d0
      cable_set=.TRUE.
      
    else
    
      do conductor1=1,conductor-2
        do conductor2=conductor1+1,conductor-1
      
          x1=xc(conductor1)
          y1=yc(conductor1)
	  r1=ri(conductor1)+delta/2d0
          x2=xc(conductor2)
          y2=yc(conductor2)
	  r2=ri(conductor2)+delta/2d0
	
	  do side=1,2
!	    write(*,*)conductor1,conductor2,side
	    call calculate_cable_centre_coordinates(x1,y1,r1,x2,y2,r2,x0,y0,r0,side)	
	  
            call check_cable_intersection (n_conductors,conductor-1,	&
	                                   xc,yc,ri,delta,		&
	                                   x0,y0,r0,		&
	                                   cable_intersection)
				       	
	    d0=sqrt(x0**2+y0**2)
	    if ( (.NOT.cable_intersection).AND.(d0.lt.dmin) ) then
	      xmin=x0
	      ymin=y0
	      dmin=d0
	      cable_set=.TRUE.
	    end if
	
	  end do ! next side
		
        end do
      end do
      
    end if ! conductor.gt.2
      
    if (.NOT.cable_set) then
      write(*,*)'Unable to set cable position'
      stop
    end if
    
    xc(conductor)=xmin
    yc(conductor)=ymin
        
    conductor=conductor+1
    
    if (conductor.le.n_conductors) goto 10
!  end do

! TLM return conductor
  conductor=n_conductors+1
  xc(conductor)=0d0
  yc(conductor)=0d0
  
! Work out the maximum extent of the bundle and shift the centre of gravity to the centre of the 
! return conductor here...

  xmin=0
  xmax=0
  ymin=0
  ymax=0
  
  do conductor=1,n_conductors
    if ( (xc(conductor)+ri(conductor)).gt.xmax) xmax=(xc(conductor)+ri(conductor))
    if ( (xc(conductor)-ri(conductor)).lt.xmin) xmin=(xc(conductor)-ri(conductor))
    if ( (yc(conductor)+ri(conductor)).gt.ymax) ymax=(yc(conductor)+ri(conductor))
    if ( (yc(conductor)-ri(conductor)).lt.ymin) ymin=(yc(conductor)-ri(conductor))
  end do
    
  xshift=(xmax+xmin)/2d0
  yshift=(ymax+ymin)/2d0
  
  do conductor=1,n_conductors
    xc(conductor)=xc(conductor)-xshift
    yc(conductor)=yc(conductor)-yshift
  end do
  
  xmax=xmax-xshift
  xmin=xmin-xshift
  ymax=ymax-yshift
  ymin=ymin-yshift
      
  rmax_found=0d0

  do conductor=1,n_conductors    
			  
    call check_bundle_dimensions( xc(conductor),yc(conductor),ri(conductor),rmax_found)
    			          
  end do ! next conductor
  
! do last (reference) conductor, set the radius (with a safety factor) to the multi-conductor cable radius
  
  conductor=n_conductors+1
    
  ri(conductor)=rmax_found*1.05d0

!  write(*,*)'FINISHED: create_bundle_cross_section_geometry'
    
  RETURN
  
  
END SUBROUTINE create_bundle_cross_section_geometry
!
! __________________________________________________
!
!  
  SUBROUTINE check_bundle_dimensions(x,y,ri,rmax_found)
!   
USE constants	
USE file_information   
USE TLM_general

IMPLICIT NONE

  real*8 x,y,ri,rmax_found
  
! local variables  

  real*8 t
  real*8 xp2,yp2,radius

  integer tloop
  
! START

! loop over theta    

  do tloop=0,50
  
    t=2d0*pi*dble(tloop)/50d0

    xp2=x+ri*cos(t)
    yp2=y+ri*sin(t)
    
    radius=sqrt((xp2*xp2)+yp2*yp2)
    
    rmax_found=max(radius,rmax_found)

  end do
  
  return
  
  END SUBROUTINE check_bundle_dimensions
!
! __________________________________________________
!
!  
  SUBROUTINE check_cable_intersection (n_conductors,last_conductor,	&
	                               xc,yc,ri,delta,		&
	                               x0,y0,r0,		&
	                               cable_intersection)

! check intersection between the proposed conductor and all
! conductors from 1 to last_condcutor

  integer n_conductors,last_conductor
  real*8 ::xc(1:n_conductors+1)
  real*8 ::yc(1:n_conductors+1)
  real*8 ::ri(1:n_conductors+1)
  real*8 delta

  real*8 x0,y0,r0
  logical cable_intersection

! local variables

  integer conductor
  real*8 separation
  real*8 min_separation

! START

  cable_intersection=.FALSE.
  
  do conductor=1,last_conductor
  
    min_separation=ri(conductor)+r0
    separation=sqrt( (x0-xc(conductor))**2+(y0-yc(conductor))**2 )
    if (separation.lt.min_separation) then
      cable_intersection=.TRUE.
      return
    end if
    
  end do

  return
  	       
  end SUBROUTINE check_cable_intersection
!
! __________________________________________________
!
!  
  SUBROUTINE calculate_cable_centre_coordinates(x1,y1,r1,x2,y2,r2,x0,y0,r0,side)
  
  real*8  x1,y1,r1,x2,y2,r2,x0,y0,r0
  integer side

! local variables

  real*8 separation
  real*8 a,b,c
  real*8 ctc
  real*8 theta,theta1,theta2

! START

  separation=sqrt( (x1-x2)**2+(y1-y2)**2 )
  
  if ( separation.ge.(r1+r2+2d0*r0) ) then
! new wire fits between the two conductors so return the midpoint
    x0=(x1+x2)/2d0
    y0=(y1+y2)/2d0
!    write(*,*)'Trial coordinates,b',x0,y0
    return
  end if
  
  a=separation
  b=r0+r1
  c=r0+r2
  
  ctc=(a*a+b*b-c*c)/(2d0*a*b)
  theta2=atan2( (y2-y1),(x2-x1) )
  theta1=acos(ctc)
  
  if (side.eq.1) then
    theta=theta2-theta1
  else
    theta=theta2+theta1
  end if
  
  x0=x1+b*cos(theta)
  y0=y1+b*sin(theta)
  
!  write(*,*)'side=',side,' ctc',ctc
!  write(*,*)theta,theta1,theta2
!  write(*,*)'Trial coordinates',x0,y0
  
  RETURN
  
  END SUBROUTINE calculate_cable_centre_coordinates
