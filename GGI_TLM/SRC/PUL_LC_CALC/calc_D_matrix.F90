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
  SUBROUTINE calc_D_matrix()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  
  
  real*8 xp,yp,rp,tp
  real*8 xw,yw,rw,tw
  real*8 xc,yc
  real*8 phi_local

  integer row,col
  integer row_point,col_point
  integer row_wire,col_wire
  integer m
  integer wire
  integer count
  
  logical TLM_return_source_flag

! START
 
! potential evaluation points go down columns,   
! Fourier terms for each wire go along rows

  count=1
  do wire=1,pul_nwires  
    pul_wire_spec(wire)%matrix_position=count
    count=count+1+2*pul_wire_spec(wire)%nterms	
  end do ! next wire
  
  pul_D(:,:)=0d0
 
  row=0
  do row_wire=1,pul_nwires 
  
    do row_point=1,pul_wire_spec(row_wire)%npoints
    
      row=row+1

! coordinates of test point      
      xp=pul_wire_spec(row_wire)%xp(row_point)
      yp=pul_wire_spec(row_wire)%yp(row_point) 
    
      col=0
      do col_wire=1,pul_nwires 
      
        if ( (col_wire.eq.pul_nwires).AND.(pul_include_TLM_return) ) then
          TLM_return_source_flag=.TRUE.
	else
	  TLM_return_source_flag=.FALSE.
	end if
	
        do col_point=1,pul_wire_spec(col_wire)%npoints
          col=col+1
	
! coordinates of wire source point
          xc=pul_wire_spec(col_wire)%xc
          yc=pul_wire_spec(col_wire)%yc
		
! calculate geometry wrt centre of the current wire

          rw= pul_wire_spec(col_wire)%rw
	
          rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
          tp=atan2( (yp-yc),(xp-xc) )
     
          if (TLM_return_source_flag) then ! the potentatial point is inside the TLM return conductor
	  
            if (col_point.eq.1) then ! constant term
	  
	      m=0
	      call calc_phi_ICR_const(rp,tp,rw,phi_local)
	      pul_D(row,col)=phi_local
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	  
	      m=col_point-1
	      call calc_phi_ICR_cos(rp,tp,rw,m,phi_local)
	      pul_D(row,col)=phi_local
	    
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	  
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_phi_ICR_sin(rp,tp,rw,m,phi_local)
	      pul_D(row,col)=phi_local          
	    
            else
	  
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	    
	    end if
     
          else ! not TLM_return_source
     
            if (col_point.eq.1) then ! constant term
	  
	      m=0
	      call calc_phi_OCR_const(rp,tp,rw,phi_local)
	      pul_D(row,col)=-rw*log(rp)/eps0
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	  
	      m=col_point-1
	      call calc_phi_OCR_cos(rp,tp,rw,m,phi_local)
	      pul_D(row,col)=phi_local
	    
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	  
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_phi_OCR_sin(rp,tp,rw,m,phi_local)
	      pul_D(row,col)=phi_local
	    
            else
	  
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	    
	    end if
	    
	  end if ! TLM_return_source_flag
	  	  
        end do ! next col point   
    
      end do ! next col wire
          
    end do ! next row point   
    
  end do ! next row wire
  
  return
  
  end SUBROUTINE calc_D_matrix
!
! ______________________________________________________
!
!  
  SUBROUTINE calc_D_matrix_dielectric()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  
  
  real*8 xp,yp,rp,tp
  real*8 xw,yw,rw,tw
  real*8 xc,yc
  real*8 rb
  
  real*8 epsr
  
  real*8 phi_local
  
  real*8 Er,Et
  real*8 Ex,Ey
  real*8 nx,ny
  real*8 En,pul_Dn_in,pul_Dn_out
  real*8 rx,ry,tx,ty

  integer row,col
  integer row_point,col_point
  integer row_wire,col_wire
  integer m
  integer wire
  integer count
  
  logical TLM_return_source_flag

! START
 
! potential evaluation points go down columns,   
! Fourier terms for each wire go along rows

  count=1
  do wire=1,pul_nwires  
    pul_wire_spec(wire)%matrix_position=count
    count=count+1+2*pul_wire_spec(wire)%nterms	
  end do ! next wire

  do wire=1,pul_nwires  
    pul_wire_spec(wire)%matrix_position2=count
    count=count+1+2*pul_wire_spec(wire)%nterms	
  end do ! next wire
  
  pul_D(:,:)=0d0

! pul_D11 SUB-MATRIX
! CONDUCTOR POTENTIAL DUE TO (FREE CHARGE ON CONDUCTORS - BOUND CHARGE ON CONDUCTORS)
 
!  write(*,*)'pul_D11 SUB-MATRIX'
 
  row=0
  do row_wire=1,pul_nwires 
  
    do row_point=1,pul_wire_spec(row_wire)%npoints
    
      row=row+1

! coordinates of test point      
      xp=pul_wire_spec(row_wire)%xp(row_point)
      yp=pul_wire_spec(row_wire)%yp(row_point) 
    
      col=0
      do col_wire=1,pul_nwires 
      
        if ( (col_wire.eq.pul_nwires).AND.(pul_include_TLM_return) ) then
          TLM_return_source_flag=.TRUE.
	else
	  TLM_return_source_flag=.FALSE.
	end if
  
        do col_point=1,pul_wire_spec(col_wire)%npoints
          col=col+1
	
! coordinates of wire source point      
          xc=pul_wire_spec(col_wire)%xc
          yc=pul_wire_spec(col_wire)%yc
			
! calculate geometry wrt centre of the current wire

          rw= pul_wire_spec(col_wire)%rw
	
          rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
          tp=atan2( (yp-yc),(xp-xc) )
	  
	  rb=pul_wire_spec(col_wire)%rw ! source charge radius is radius of conductor
	  
          if (TLM_return_source_flag) then ! THE POTENTATIAL POINT IS INSIDE THE TLM RETURN CONDUCTOR USE FORMULA TYPE B
	  
            if (col_point.eq.1) then ! constant term
	  
	      m=0
	      call calc_phi_ICR_const(rp,tp,rb,phi_local)
	      pul_D(row,col)=phi_local
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	  
	      m=col_point-1
	      call calc_phi_ICR_cos(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local
	    
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	  
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_phi_ICR_sin(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local          
	    
            else
	  
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	    
	    end if
     
          else ! not TLM_return_source
! EVALUATION POINT IS OUTSIDE CHARGE DENSITY SO USE FORMULA TYPE A
     
            if (col_point.eq.1) then ! constant term
	  
	      m=0
	      call calc_phi_OCR_const(rp,tp,rb,phi_local)
	      pul_D(row,col)=phi_local
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	  
	      m=col_point-1
	      call calc_phi_OCR_cos(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local
	  
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	  
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_phi_OCR_sin(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local
	  
            else
	  
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	  
	    end if
	    
	  end if ! TLM_return_source
	    	  	  
        end do ! next col point   
    
      end do ! next col wire
          
    end do ! next row point   
    
  end do ! next row wire
 
! pul_D12 SUB-MATRIX
! CONDUCTOR POTENTIAL DUE TO BOUND CHARGE ON DIELECTRICS
 
!  write(*,*)'pul_D12 SUB-MATRIX'
 
  row=0
  do row_wire=1,pul_nwires
  
    do row_point=1,pul_wire_spec(row_wire)%npoints
    
      row=row+1

! coordinates of test point      
      xp=pul_wire_spec(row_wire)%xp(row_point)
      yp=pul_wire_spec(row_wire)%yp(row_point) 
    
      col=pul_tot_npoints
      do col_wire=1,pul_nwires_in  ! note, don't include dielectric terms for TLM return
  
        do col_point=1,pul_wire_spec(col_wire)%npoints
          col=col+1
	
! coordinates of wire source point      
          xc=pul_wire_spec(col_wire)%xc
          yc=pul_wire_spec(col_wire)%yc
		
! calculate geometry wrt centre of the current wire

          rw= pul_wire_spec(col_wire)%rw
	
          rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
          tp=atan2( (yp-yc),(xp-xc) )
	  
	  rb=pul_wire_spec(col_wire)%ri ! source charge radius is radius of dielectric   
	  
	  if (row_wire.ne.col_wire) then ! EVALUATION POINT IS OUTSIDE CHARGE DENSITY SO USE FORMULA TYPE A
     
            if (col_point.eq.1) then ! constant term
	  
	      m=0
	      call calc_phi_OCR_const(rp,tp,rb,phi_local)
	      pul_D(row,col)=phi_local
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	  
	      m=col_point-1
	      call calc_phi_OCR_cos(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local
	    
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	  
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_phi_OCR_sin(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local
	    
            else
	  
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	    
	    end if
	    
	  else ! EVALUATION POINT IS INSIpul_DE CHARGE DENSITY SO USE FORMULA TYPE B
	  
            if (col_point.eq.1) then ! constant term
	  
	      m=0
	      call calc_phi_ICR_const(rp,tp,rb,phi_local)
	      pul_D(row,col)=phi_local
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	  
	      m=col_point-1
	      call calc_phi_ICR_cos(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local
	    
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	  
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_phi_ICR_sin(rp,tp,rb,m,phi_local)
	      pul_D(row,col)=phi_local          
	    
            else
	  
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	    
	    end if
	    
	  end if ! source point inside bound charge region
	  	  
        end do ! next col point   
    
      end do ! next col wire
          
    end do ! next row point   
    
  end do ! next row wire
  
! pul_D21 SUB-MATRIX
! NORMAL D FIELD ON DIELECTRIC DUE TO CHARGE ON CONDUCTORS
 
!  write(*,*)'pul_D21 SUB-MATRIX'
 
  row=pul_tot_npoints
  do row_wire=1,pul_nwires_in  ! note, don't include dielectric terms for TLM return 
  
    do row_point=1,pul_wire_spec(row_wire)%npoints
    
      row=row+1

! coordinates of test point      
      xp=pul_wire_spec(row_wire)%xdp(row_point)
      yp=pul_wire_spec(row_wire)%ydp(row_point) 
      
! normal to dielectric at test point      
      nx=pul_wire_spec(row_wire)%nxdp(row_point)
      ny=pul_wire_spec(row_wire)%nydp(row_point)
      
      epsr=pul_wire_spec(row_wire)%epsr
    
      col=0
      do col_wire=1,pul_nwires 
      
        if ( (col_wire.eq.pul_nwires).AND.(pul_include_TLM_return) ) then
          TLM_return_source_flag=.TRUE.
	else
	  TLM_return_source_flag=.FALSE.
	end if
  
        do col_point=1,pul_wire_spec(col_wire)%npoints
          col=col+1
	
! coordinates of wire source point      
          xc=pul_wire_spec(col_wire)%xc
          yc=pul_wire_spec(col_wire)%yc
		
! calculate geometry wrt centre of the current wire
	
          rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
          tp=atan2( (yp-yc),(xp-xc) )
	  
! calculate projections of r and theta directions onto the x and y directions
          rx= cos(tp)
	  ry= sin(tp)
	  tx=-sin(tp)
	  ty= cos(tp)
	  
	  rb=pul_wire_spec(col_wire)%rw ! source charge radius is radius of conductor   
	  
          if (TLM_return_source_flag) then ! THE FIELD POINT IS INSIDE THE TLM RETURN CONDUCTOR USE FORMULA TYPE B
     
            if (col_point.eq.1) then ! constant term
	 
	      m=0
	      call calc_E_ICR_const(rp,tp,rb,Er,Et)
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	 
	      m=col_point-1
	      call calc_E_ICR_cos(rp,tp,rb,m,Er,Et)
	  
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	 
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_E_ICR_sin(rp,tp,rb,m,Er,Et)
	  
            else
	 
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	  
	    end if

          else ! EVALUATION POINT OUTSIDE CHARGE DENSITY SO USE FORMULA TYPE A
     
            if (col_point.eq.1) then ! constant term
	 
	      m=0
	      call calc_E_OCR_const(rp,tp,rb,Er,Et)
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	 
	      m=col_point-1
	      call calc_E_OCR_cos(rp,tp,rb,m,Er,Et)
	  
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	 
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_E_OCR_sin(rp,tp,rb,m,Er,Et)
	  
            else
	 
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	  
	    end if
	  
	  end if

! calculate E field in cartesian coordinates
          Ex=Er*rx+Et*tx
	  Ey=Er*ry+Et*ty
	  
! project E field onto normal to dielectric
	  En=Ex*nx+Ey*ny
	  pul_Dn_in=En*epsr
	  pul_Dn_out=En

! jump in pul_D due to charges on conductor surfaces only	  
          
	  pul_D(row,col)=pul_Dn_in-pul_Dn_out
	  	  
        end do ! next col point   
    
      end do ! next col wire
          
    end do ! next row point   
    
  end do ! next row wire
   
! pul_D22 SUB-MATRIX
! NORMAL pul_D FIELpul_D ON pul_DIELECTRIC pul_DUE TO BOUNpul_D CHARGE ON pul_DIELECTRICS
 
!  write(*,*)'pul_D22 SUB-MATRIX'
 
  row=pul_tot_npoints
  do row_wire=1,pul_nwires_in  ! note, don't include dielectric terms for TLM return 
  
    do row_point=1,pul_wire_spec(row_wire)%npoints
    
      row=row+1

! coordinates of test point      
      xp=pul_wire_spec(row_wire)%xdp(row_point)
      yp=pul_wire_spec(row_wire)%ydp(row_point) 
      
! normal to dielectric at test point      
      nx=pul_wire_spec(row_wire)%nxdp(row_point)
      ny=pul_wire_spec(row_wire)%nydp(row_point)
      
      epsr=pul_wire_spec(row_wire)%epsr
    
      col=pul_tot_npoints
      do col_wire=1,pul_nwires_in  ! note, don't include dielectric terms for TLM return 
  
        do col_point=1,pul_wire_spec(col_wire)%npoints
          col=col+1
	
! coordinates of wire source point      
          xc=pul_wire_spec(col_wire)%xc
          yc=pul_wire_spec(col_wire)%yc
		
! calculate geometry wrt centre of the current wire

          rw= pul_wire_spec(col_wire)%rw
	
          rp= sqrt( (xp-xc)**2+(yp-yc)**2 )
          tp=atan2( (yp-yc),(xp-xc) )
	  
! calculate projections of r and theta directions onto the x and y directions
          rx= cos(tp)
	  ry= sin(tp)
	  tx=-sin(tp)
	  ty= cos(tp)
	  
	  rb=pul_wire_spec(col_wire)%ri ! source charge radius is radius of dielectric  
	  
	  if (row_wire.ne.col_wire) then ! EVALUATION POINT IS OUTSIDE CHARGE DENSITY SO USE FORMULA TYPE A FOR IN AND OUT
     
            if (col_point.eq.1) then ! constant term
	 
	      m=0
	      call calc_E_OCR_const(rp,tp,rb,Er,Et)
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	 
	      m=col_point-1
	      call calc_E_OCR_cos(rp,tp,rb,m,Er,Et)
	  
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	 
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_E_OCR_sin(rp,tp,rb,m,Er,Et)
	  
            else
	 
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	  
	    end if

! calculate E field in cartesian coordinates
            Ex=Er*rx+Et*tx
	    Ey=Er*ry+Et*ty
	  
! project E field onto normal to dielectric
	    En=Ex*nx+Ey*ny
	    pul_Dn_in=En*epsr
	    pul_Dn_out=En

! jump in pul_D due to charges on conductor surfaces only	  
	    pul_D(row,col)=pul_Dn_in-pul_Dn_out
	    
	  else ! EVALUATION POINT LIES ON THIS DIELECTRIC INTERFACE SO USE DIFFERENT FORMLAE FOR INSIDE AND OUTSIDE

! INSIDE EVALUATION POINT IS INSIDE THIS CHARGE DENSITY SO USE FORMULA TYPE B 
	  
            if (col_point.eq.1) then ! constant term
	 
	      m=0
	      call calc_E_ICR_const(rp,tp,rb,Er,Et)
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	 
	      m=col_point-1
	      call calc_E_ICR_cos(rp,tp,rb,m,Er,Et)
	  
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	 
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_E_ICR_sin(rp,tp,rb,m,Er,Et)
	  
            else
	 
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	  
	    end if

! calculate E field in cartesian coordinates
            Ex=Er*rx+Et*tx
	    Ey=Er*ry+Et*ty
	  
! project E field onto normal to dielectric and calculate pul_D_in
	    En=Ex*nx+Ey*ny
	    pul_Dn_in=En*epsr
	    
! OUTSIDE EVALUATION POINT IS OUTSIDE CHARGE DENSITY SO USE FORMULA TYPE A FOR OUTSIDE
     
            if (col_point.eq.1) then ! constant term
	 
	      m=0
	      call calc_E_OCR_const(rp,tp,rb,Er,Et)
    
	    else if (col_point.le.1+pul_wire_spec(col_wire)%nterms) then ! A term (cos theta)
	 
	      m=col_point-1
	      call calc_E_OCR_cos(rp,tp,rb,m,Er,Et)
	  
	    else if (col_point.le.1+2*pul_wire_spec(col_wire)%nterms) then ! B term (sin theta)
	 
	      m=col_point-pul_wire_spec(col_wire)%nterms-1
	      call calc_E_OCR_sin(rp,tp,rb,m,Er,Et)
	  
            else
	 
	      write(*,*)'Matrix fill error'
	      write(*,*)'col_wire=',col_wire
	      write(*,*)'col_point=',col_point
	      write(*,*)'row_wire=',row_wire
	      write(*,*)'row_point=',row_point
	      stop
	  
	    end if

! calculate E field in cartesian coordinates
            Ex=Er*rx+Et*tx
	    Ey=Er*ry+Et*ty
	  
! project E field onto normal to dielectric and calculate pul_D_out
	    En=Ex*nx+Ey*ny
	    pul_Dn_out=En
	 
! jump in pul_D due to charges on conductor surfaces only	  
	    pul_D(row,col)=pul_Dn_in-pul_Dn_out	  
	  
	  end if ! point inside charge density volume	 
	   
        end do ! next col point   
    
      end do ! next col wire
          
    end do ! next row point   
    
  end do ! next row wire

  return
  
  end SUBROUTINE calc_D_matrix_dielectric
  
