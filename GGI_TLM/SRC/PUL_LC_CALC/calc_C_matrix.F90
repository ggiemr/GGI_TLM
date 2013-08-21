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
  SUBROUTINE calc_C_matrix()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  
  
  integer row,col
  integer row_wire,col_wire
  integer Crow,Ccol
  
  real*8 sum,row_sum,col_sum,tot_sum

! START
  
! calculate Capacitance matrix
!  write(*,*)'Calculate C'

  tot_sum=0d0
  do row=1,pul_nwires
    do col=1,pul_nwires
      tot_sum=tot_sum+pul_Cg(row,col)
    end do
  end do

  Crow=0
  do row_wire=1,pul_nwires
  
    if (row_wire.ne.pul_return_conductor) then
      Crow=Crow+1
    
      Ccol=0
      do col_wire=1,pul_nwires
      
        if (col_wire.ne.pul_return_conductor) then
        
	  Ccol=Ccol+1
	
          row_sum=0d0
          do row=1,pul_nwires
            row_sum=row_sum+pul_Cg(row,col_wire)
          end do
      
          col_sum=0d0
          do col=1,pul_nwires
            col_sum=col_sum+pul_Cg(row_wire,col)
          end do
            
          pul_C(Crow,Ccol)=pul_Cg(row_wire,col_wire)-(row_sum*col_sum)/tot_sum
	  
        end if ! col not return conductor
	
      end do ! next col_wire
      
    end if ! row not return conductor
    
  end do ! next row_wire
  
  return
  
  end SUBROUTINE calc_C_matrix
  
