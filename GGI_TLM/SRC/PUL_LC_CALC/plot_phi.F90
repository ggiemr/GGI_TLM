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
  SUBROUTINE plot_phi()
   
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
  
  real*8 phi_local,phi_temp
  real*8 Ex,Ey
  real*8 Dxp,Dyp
  real*8 epsr
  
  real*8,allocatable :: phi_save(:,:)

! START

  dx=(pul_xmax-pul_xmin)/(pul_nxp-1)
  
  dy=(pul_ymax-pul_ymin)/(pul_nyp-1)
  
  allocate ( phi_save(1:pul_nxp,1:pul_nyp) )

! 2d potential plot 
 
  open(unit=20,file='phi2d_1.dat')
  
  do xloop=1,pul_nxp
  
    xp=pul_xmin+(xloop-1)*dx
  
    do yloop=1,pul_nyp
    
      yp=pul_ymin+(yloop-1)*dy
      
      phi_local=0d0

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

! loop over terms relating to this wire	
        do term=1,pul_wire_spec(wire)%npoints
	
	  row=row+1
     
          if (term.eq.1) then ! constant term
	  
            m=0
	    call calc_phi_const(rp,tp,rw,phi_temp)
	    phi_local=phi_local+phi_temp*pul_alpha(row)
	    
	  else if (term.le.1+pul_wire_spec(wire)%nterms) then ! A term (cos theta)
	  
	    m=term-1
	    call calc_phi_cos(rp,tp,rw,m,phi_temp)
	    phi_local=phi_local+phi_temp*pul_alpha(row)
	    	    
	  else if (term.le.1+2*pul_wire_spec(wire)%nterms) then ! B term (sin theta)
	  
	    m=term-pul_wire_spec(wire)%nterms-1
	    call calc_phi_sin(rp,tp,rw,m,phi_temp)
	    phi_local=phi_local+phi_temp*pul_alpha(row)	    
	    	    
	  end if ! term type
	    
	end do ! next term
     
      end do ! next wire
      
      if (pul_include_dielectric) then
      
! loop over wires      
     
      do wire=1,pul_nwires_in
      		
! coordinates of wire source point      
        xc=pul_wire_spec(wire)%xc
        yc=pul_wire_spec(wire)%yc
		
! calculate geometry wrt centre of the current wire

        rw= pul_wire_spec(wire)%ri  ! radius of insulator 
	
        rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
        tp=atan2( (yp-yc),(xp-xc) )

! loop over terms relating to this wire	
        do term=1,pul_wire_spec(wire)%npoints
	
	  row=row+1
     
          if (term.eq.1) then ! constant term
	  
            m=0
	    call calc_phi_const(rp,tp,rw,phi_temp)
	    phi_local=phi_local+phi_temp*pul_alpha(row)
	    
	  else if (term.le.1+pul_wire_spec(wire)%nterms) then ! A term (cos theta)
	  
	    m=term-1
	    call calc_phi_cos(rp,tp,rw,m,phi_temp)
	    phi_local=phi_local+phi_temp*pul_alpha(row)
	    	    
	  else if (term.le.1+2*pul_wire_spec(wire)%nterms) then ! B term (sin theta)
	  
	    m=term-pul_wire_spec(wire)%nterms-1
	    call calc_phi_sin(rp,tp,rw,m,phi_temp)
	    phi_local=phi_local+phi_temp*pul_alpha(row)	    
	    	    
	  end if ! term type
	    
	end do ! next term
     
      end do ! next wire     
      
      end if
      
      write(20,8000)xp,yp,phi_local
8000  format(3E16.7)

      phi_save(xloop,yloop)=phi_local
    
    end do ! next y
    
    write(20,*)
    
  end do ! next x
  
  close(unit=20)
 
 
! 2d E  field plot 
 
  open(unit=20,file='E2d_2.dat')
  open(unit=21,file='D2d_2.dat')
  open(unit=22,file='epsr.dat')
  
  do xloop=2,pul_nxp-1
  
    xp=pul_xmin+(xloop-1)*dx
  
    do yloop=2,pul_nyp-1
    
      yp=pul_ymin+(yloop-1)*dy
      
      epsr=1d0
! loop over wires      
      do wire=1,pul_nwires_in
      		
! coordinates of wire source point      
        xc=pul_wire_spec(wire)%xc
        yc=pul_wire_spec(wire)%yc
		
! calculate geometry wrt centre of the current wire

        rw= pul_wire_spec(wire)%ri
	
        rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
	  
	if (rp.lt.rw) epsr=pul_wire_spec(wire)%epsr	
      
      end do ! next wire
      
      Ex=-(phi_save(xloop+1,yloop)-phi_save(xloop-1,yloop))/(2d0*dx)
      Ey=-(phi_save(xloop,yloop+1)-phi_save(xloop,yloop-1))/(2d0*dy)
      
      Dxp=Ex*epsr
      Dyp=Ey*epsr
      
      write(20,8010)xp,yp,Ex,Ey
      write(21,8010)xp,yp,Dxp,Dyp
      write(22,8010)xp,yp,epsr
      
8010  format(4E16.7)
    
    end do ! next y
    
    write(20,*)
    write(21,*)
    write(22,*)
    
  end do ! next x
  
  close(unit=20)
  close(unit=21)
  close(unit=22)
 
  deallocate ( phi_save )
 
 return
 
 END SUBROUTINE plot_phi
