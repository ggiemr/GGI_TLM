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
!       subroutine cmatmul(A,ar,ac,B,br,bc,C,matdim)
!       subroutine ctranspose(u,nr,nc,ut,nmodes)
!       subroutine cmatvmul(A,ar,ac,B,br,C,matdim)
!       subroutine csvd(A,ar,ac,matdim,U,GAMMA,VT) 
!       subroutine csvd_invert(A,ar,ac,AI,matdim) 
!
! __________________________________________________
!
!
       subroutine ctranspose(u,nr,nc,ut,nmodes)

IMPLICIT NONE
       
       integer nmodes
       complex*16 u(nmodes,nmodes),ut(nmodes,nmodes)
       integer nr,nc
       integer r,c
       
       do r=1,nr
         do c=1,nc
	   ut(c,r)=conjg(u(r,c))
	 end do
       end do

       return
       end
!
! __________________________________________________
!
!
       subroutine  cmatmul(A,ar,ac,B,br,bc,C,matdim)

IMPLICIT NONE

       
       integer ar,ac,br,bc,matdim
       complex*16 A(matdim,matdim)
       complex*16 B(matdim,matdim)
       complex*16 C(matdim,matdim)
       
       complex*16 sum
       
       integer row,col,i
  
       if (ac.ne.br) then
         print*,'matrix dimension error in cmatmul'
	 print*,'ac=',ac
	 print*,'br=',br
	 stop
       end if
  
       do row=1,ar
         do col=1,bc
	   c(row,col)=0.0
           do i=1,ac
	     c(row,col)=c(row,col)+a(row,i)*b(i,col)
	   end do
	 end do
       end do
             
       return
       end
!
! __________________________________________________
!
!
       subroutine  cmatvmul(A,ar,ac,B,br,C,matdim)

IMPLICIT NONE

! multiply complex vector by a complex matrix
       
       integer ar,ac,br,matdim
       complex*16 A(matdim,matdim)
       complex*16 B(matdim)
       complex*16 C(matdim)
       
       integer i,j
       complex*16 sum
              
       if (br.ne.ac) then
         print*,'dimension error in cmatvmul'
	 print*,'number of columns in a=',ac
	 print*,'number of rows in b   =',br
	 stop
       end if
         
       do i=1,ar
         sum=(0.0,0.0)
         do j=1,ac
	   sum=sum+A(i,j)*B(j)
	 end do
	 C(i)=sum
       end do	 
       
       return
       end
!
! __________________________________________________
!
! 
  SUBROUTINE cinvert_Gauss_Jordan(A,n,AI,matdim) 

IMPLICIT NONE

! invert matrix using Gauss Jordan method with pivoting
       
  integer	:: n,matdim
  complex*16	:: A(matdim,matdim)
  complex*16	:: AI(matdim,matdim)

! local variables
  		       
  integer	:: row,col,reduce_col,i
  
  real*8	:: max_element
  complex*16	:: pivot_element
  integer	:: max_row
  
  integer	:: pivot_row
  
  integer	:: pivot_row_save(matdim)
  
  complex*16	:: row_multiplier
  complex*16	:: swap

! START
  
! copy A to AI 
  AI(1:n,1:n)= A(1:n,1:n)

  pivot_row_save(1:matdim)=0
  
! loop over columns of the matrix and reduce each column in turn to identity matrix column     
  do reduce_col=1,n
  
! find the largest element in this column and use as the pivot element
    max_element=0d0
    max_row=0
    do row=reduce_col,n
      if (abs(AI(row,reduce_col)).GT.max_element) then
        max_element=abs(AI(row,reduce_col))
	max_row=row
      end if
    end do  
    
    if (max_row.eq.0) then
! all elements are zero so singular matrix
      write(*,*)'Singular matrix found in cinvert_Gauss_Jordan'
      STOP
    end if
    
    pivot_row=max_row
    pivot_row_save(reduce_col)=pivot_row
    
! swap pivot row with the row reduce_col

    if (pivot_row.ne.reduce_col) then
      do col=1,n
        swap=AI(reduce_col,col)
	AI(reduce_col,col)=AI(pivot_row,col)
	AI(pivot_row,col)=swap
      end do 
    end if
    
    pivot_row=reduce_col   
    pivot_element=AI(reduce_col,reduce_col)
    
! operate on pivot row    
    do col=1,n
      if (col.ne.reduce_col) then
        AI(pivot_row,col) = AI(pivot_row,col)/pivot_element
      else	
        AI(pivot_row,col) = (1d0,0d0)/pivot_element
      end if
    end do

! operate on rows other than the pivot row   
    do row=1,n
    
      if (row.ne.pivot_row) then
      
        row_multiplier=AI(row,reduce_col)
	
        do col=1,n
          if (col.ne.reduce_col) then
            AI(row,col) = AI(row,col)- AI(pivot_row,col)*row_multiplier
          else	
            AI(row,reduce_col) =-AI(pivot_row,reduce_col)*row_multiplier
          end if
        end do
	
      end if ! not pivot row
    
    end do ! next row
    
  end do ! next column of the matrix to reduce
  
  do reduce_col=n,1,-1
  
    if (reduce_col.ne.pivot_row_save(reduce_col)) then
! rows were swapped so must swap the corresponding columns

      do row=1,n
        swap=AI(row,pivot_row_save(reduce_col))
        AI(row,pivot_row_save(reduce_col))=AI(row,reduce_col)
	AI(row,reduce_col)=swap
      end do
      
    end if
    
  end do

  RETURN
  
  END
!
! __________________________________________________
!
! 
  SUBROUTINE cinvert_Moore_Penrose(A,m,n,AI,matdim) 

IMPLICIT NONE

! invert matrix using the Morse Penrose generalised inverse - invert with Gauss Jordan method
       
  integer	:: m,n,matdim
  complex*16	:: A(matdim,matdim)
  complex*16	:: AI(matdim,matdim)

! local variables
  complex*16	:: AH(matdim,matdim)
  complex*16	:: T(matdim,matdim)
  complex*16	:: TI(matdim,matdim)
  		       
! START

  if (m.ge.n) then
  
    CALL ctranspose(A,m,n,AH,matdim)
    CALL cmatmul(AH,n,m,A,m,n,T,matdim)
    CALL cinvert_Gauss_Jordan(T,n,TI,matdim) 
    CALL cmatmul(TI,n,n,AH,n,m,AI,matdim)
    
  else
  
    CALL ctranspose(A,m,n,AH,matdim)
    CALL cmatmul(A,m,n,AH,n,m,T,matdim)
    CALL cinvert_Gauss_Jordan(T,m,TI,matdim) 
    CALL cmatmul(AH,n,m,TI,m,m,AI,matdim)
  
  end if
  

  RETURN
  
  END
