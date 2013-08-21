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
  SUBROUTINE plot_E()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  
  
  real*8 dx
  real*8 dy
  integer xloop,yloop
  
  integer wire,term
  integer row
  integer m
  
  real*8 xp,yp,rp,tp
  real*8 xw,yw,rw
  real*8 xc,yc
  
  real*8 Ex_local,Ey_local
  real*8 Er_temp,Et_temp
 
  real*8 Ex,Ey
  real*8 epsEx,epsEy
  real*8 epsr
  real*8 conductor_constant
  real*8 rx,ry,tx,ty

! START

  dx=(pul_xmax-pul_xmin)/(pul_nxp-1)
  
  dy=(pul_ymax-pul_ymin)/(pul_nyp-1)

! 2d plot 
 
  open(unit=20,file='E2d_1.dat')
  open(unit=21,file='D2d_1.dat')
  
  do xloop=1,pul_nxp
  
    xp=pul_xmin+(xloop-1)*dx
  
    do yloop=1,pul_nyp
    
      yp=pul_ymin+(yloop-1)*dy
      
      epsr=1d0
      conductor_constant=1d0
      
      Ex_local=0d0
      Ey_local=0d0

! loop over wires      
      row=0
      do wire=1,pul_nwires 
      		
! coordinates of wire source point      
        xc=pul_wire_spec(wire)%xc
        yc=pul_wire_spec(wire)%yc
		
! calculate geometry wrt centre of the current wire

        rw= pul_wire_spec(wire)%rw
	
        rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
        tp=atan2( (yp-yc),(xp-xc) )
	
	if ( (wire.ne.pul_nwires).AND.(rp.lt.rw) ) then
	  conductor_constant=0d0	
	end if
	if ( (wire.eq.pul_nwires).AND.(.NOT.pul_include_TLM_return).AND.(rp.lt.rw) ) then
	  conductor_constant=0d0	
	end if
	
! calculate projections of r and theta directions onto the x and y directions
        rx= cos(tp)
	ry= sin(tp)
	tx=-sin(tp)
	ty= cos(tp)

! loop over terms relating to this wire	
        do term=1,pul_wire_spec(wire)%npoints
	
	  row=row+1
     
          if (term.eq.1) then ! constant term
	  
            m=0
	    call calc_E_const(rp,tp,rw,Er_temp,Et_temp)
	    
	  else if (term.le.1+pul_wire_spec(wire)%nterms) then ! A term (cos theta)
	  
	    m=term-1
	    call calc_E_cos(rp,tp,rw,m,Er_temp,Et_temp)
	    	    
	  else if (term.le.1+2*pul_wire_spec(wire)%nterms) then ! B term (sin theta)
	  
	    m=term-pul_wire_spec(wire)%nterms-1
	    call calc_E_sin(rp,tp,rw,m,Er_temp,Et_temp)
	    	    
	  end if ! term type
	  
	  Ex_local=Ex_local+(Er_temp*rx+Et_temp*tx)*pul_alpha(row)
	  Ey_local=Ey_local+(Er_temp*ry+Et_temp*ty)*pul_alpha(row)
	    
	end do ! next term
     
      end do ! next wire
           
      if (pul_include_dielectric) then

! loop over wires      
        do wire=1,pul_nwires_in
      		
! coordinates of wire source point      
          xc=pul_wire_spec(wire)%xc
          yc=pul_wire_spec(wire)%yc
		
! calculate geometry wrt centre of the current wire

          rw= pul_wire_spec(wire)%ri
	
          rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
          tp=atan2( (yp-yc),(xp-xc) )
	  
	  if (rp.lt.rw) epsr=pul_wire_spec(wire)%epsr	
	
! calculate projections of r and theta directions onto the x and y directions
          rx= cos(tp)
	  ry= sin(tp)
	  tx=-sin(tp)
	  ty= cos(tp)

! loop over terms relating to this wire	
          do term=1,pul_wire_spec(wire)%npoints
	
	    row=row+1
     
            if (term.eq.1) then ! constant term
	  
              m=0
	      call calc_E_const(rp,tp,rw,Er_temp,Et_temp)
	    
	    else if (term.le.1+pul_wire_spec(wire)%nterms) then ! A term (cos theta)
	  
	      m=term-1
	      call calc_E_cos(rp,tp,rw,m,Er_temp,Et_temp)
	    	    
	    else if (term.le.1+2*pul_wire_spec(wire)%nterms) then ! B term (sin theta)
	  
	      m=term-pul_wire_spec(wire)%nterms-1
	      call calc_E_sin(rp,tp,rw,m,Er_temp,Et_temp)
	    	    
	    end if ! term type
	  
	    Ex_local=Ex_local+(Er_temp*rx+Et_temp*tx)*pul_alpha(row)
	    Ey_local=Ey_local+(Er_temp*ry+Et_temp*ty)*pul_alpha(row)
	  
	  end do ! next term
     
        end do ! next wire    
      
      end if ! dielectric  
         
! calculate E field in cartesian coordinates
      Ex=Ex_local
      Ey=Ey_local
      
      Ex=Ex*conductor_constant
      Ey=Ey*conductor_constant
      epsEx=Ex*epsr
      epsEy=Ey*epsr
	  
      write(20,8000)xp,yp,Ex,Ey
      write(21,8000)xp,yp,epsEx,epsEy
8000  format(4E14.5)
    
    end do ! next y
    
    write(20,*)
    write(21,*)
    
  end do ! next x
  
  close(unit=20)
  close(unit=21)
 
 
 return
 
 END SUBROUTINE plot_E
