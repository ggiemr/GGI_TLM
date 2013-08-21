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
  SUBROUTINE calc_Cg_matrix()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  
  
  integer row,col
  integer row_wire,col_wire
  integer row_point,col_point
  
  real*8 sum,row_sum,col_sum,tot_sum

! START
 
! extract elements of generalised capacitance matrix
!  write(*,*)'Calculate Cg'
 
  do row_wire=1,pul_nwires 
  
    row=pul_wire_spec(row_wire)%matrix_position
    
    col=0
    do col_wire=1,pul_nwires 
  
      sum=0d0
      do col_point=1,pul_wire_spec(col_wire)%npoints
        col=col+1
	     
	sum=sum+pul_B(row,col)
	    	    	  
      end do ! next col point   
      
      pul_Cg(row_wire,col_wire)=sum*2d0*pi*pul_wire_spec(row_wire)%rw

    end do ! next col wire
          
  end do ! next row wire
  
  if (pul_include_dielectric) then ! we must include bound charge terms
  
    do row_wire=1,pul_nwires_in 
  
      row=pul_wire_spec(row_wire)%matrix_position2
    
      col=0
      do col_wire=1,pul_nwires 
  
        sum=0d0
        do col_point=1,pul_wire_spec(col_wire)%npoints
          col=col+1
	     
	  sum=sum+pul_B(row,col)
	    	    	  
        end do ! next col point   
      
        pul_Cg(row_wire,col_wire)=pul_Cg(row_wire,col_wire)+sum*2d0*pi*pul_wire_spec(row_wire)%ri  ! use dielectric radius here

      end do ! next col wire
          
    end do ! next row wire
  
  end if

  return
  
  end SUBROUTINE calc_Cg_matrix
  
