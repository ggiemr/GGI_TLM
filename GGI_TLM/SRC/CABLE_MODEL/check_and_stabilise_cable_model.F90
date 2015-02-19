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
!       subroutine LC_checks
!       subroutine LC_propagation_checks
!       subroutine TLM_Cable_checks 
      
!
! Name LC_checks
!     
!
! Description
!     checks on MTL L and C matrices
!
! Comments:
!      
!
! History
!
!     started 17/03/09 CJS
!     4/02/15 CJS  Calculate the per-unit-length L and C matrices and do checks on these. 
! 		`  Write the mode velocities and impedances to file. 
!		   Also calculate the L and C matrices for the internal conductor system i.e. removing the 
!		   TLM reference conductor (C.Paul, Analysis of multiconductor transmission lines" p 
!			  
!

       SUBROUTINE LC_checks(L_in,C_in,length,nwires)

USE file_information
USE TLM_general
USE mesh

IMPLICIT NONE
       
       integer nwires
       real*8 L_in(nwires,nwires)
       real*8 C_in(nwires,nwires)
       real*8 length
       
       real*8 L(nwires,nwires)
       real*8 C(nwires,nwires)
       real*8 C01(nwires,nwires)
       
       real*8 L2(nwires-1,nwires-1)
       real*8 C02(nwires-1,nwires-1)
       real*8 C2(nwires-1,nwires-1)
       
       real*8 sum_row_C01,sum_col_C01,sum_total_C01
       real*8 sum_row_C,sum_col_C,sum_total_C
       
       real*8 LI(nwires,nwires)
       real*8 CI(nwires,nwires)
       integer i,j
       integer row,col
       logical symmetric_matrix
       
! START       
!       CALL write_line('CALLED: LC_checks',0,output_to_screen_flag)

! Copy the input matrices to local matrices
       L(:,:)=L_in(:,:)
       C(:,:)=C_in(:,:)

       write(cable_info_file_unit,*)'C:'
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(C(i,j),j=1,nwires)
       end do
       call dsvd_invert(C,nwires,nwires,Ci,nwires) 
       write(cable_info_file_unit,*)'Ci:'
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(Ci(i,j),j=1,nwires)
       end do
       
       write(cable_info_file_unit,*)'L:'
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(L(i,j),j=1,nwires)
       end do
       call dsvd_invert(L,nwires,nwires,Li,nwires) 
       write(cable_info_file_unit,*)'Li:'
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(Li(i,j),j=1,nwires)
       end do
       
       write(cable_info_file_unit,*)' '
       write(cable_info_file_unit,*)'LC checks'
       
       write(cable_info_file_unit,*)'Nwires=',nwires,' nwires=',nwires

       write(cable_info_file_unit,*)'C='
       write(cable_info_file_unit,*)'Check symmetry C'
       call dcheck_symmetry(C,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag:',symmetric_matrix
       
       write(cable_info_file_unit,*)'Ci='
       write(cable_info_file_unit,*)'Check symmetry Ci'
       call dcheck_symmetry(Ci,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag:',symmetric_matrix
       
       write(cable_info_file_unit,*)'L='
       write(cable_info_file_unit,*)'Check symmetry L'
       call dcheck_symmetry(L,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag:',symmetric_matrix
       
       write(cable_info_file_unit,*)'Li='
       write(cable_info_file_unit,*)'Check symmetry Li'
       call dcheck_symmetry(Li,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag:',symmetric_matrix

8000   format (I4,100E14.6)
       
       write(cable_info_file_unit,*)'________________________________'
       write(cable_info_file_unit,*)''
       write(cable_info_file_unit,*)'Checks on TLM per-unit-length matrices'
       write(cable_info_file_unit,*)'________________________________'

! Copy the input matrices to local matrices
       L(:,:)=L_in(:,:)/length
       C(:,:)=C_in(:,:)/length
       
       CALL LC_propagation_checks(L,C,nwires)
       
       write(cable_info_file_unit,*)'________________________________'
       write(cable_info_file_unit,*)''
       write(cable_info_file_unit,*)'Remove reference conductor and do'
       write(cable_info_file_unit,*)'Checks on reduced TLM per-unit-length matrices'
       write(cable_info_file_unit,*)'________________________________'

       call dsvd_invert(L,nwires,nwires,Li,nwires) 
       C01(:,:)=Li(:,:)*8.988D16

! Reduce the capacitance matrices to eliminate the reference conductor 
! i.e. only allow differential modes in the TLM cell
       
       sum_total_C=0d0
       sum_total_C01=0d0
       
       do row=1,nwires
         do col=1,nwires
           sum_total_C  =sum_total_C+  C(row,col)
           sum_total_C01=sum_total_C01+C01(row,col)
         end do
       end do

       do row=1,nwires-1
         do col=1,nwires-1
	 
           sum_row_C01=0d0
	   sum_col_C01=0d0
	   sum_row_C=0d0
	   sum_col_C=0d0
	   
	   do i=1,nwires
	   
             sum_row_C01=sum_row_C01+C01(row,i)
	     sum_col_C01=sum_col_C01+C01(i,col)
	     sum_row_C  =sum_row_C  +C(row,i) 
	     sum_col_C  =sum_col_C  +C(i,col) 
	     
	   end do
	   
	   C02(row,col)=C01(row,col)-sum_row_C01*sum_col_C01/sum_total_C01
	   C2(row,col) =C(row,col)  -sum_row_C  *sum_col_C  /sum_total_C
	   
         end do
       end do

       call dsvd_invert(C02,nwires-1,nwires-1,L2,nwires-1) 
       L2(:,:)=L2(:,:)*8.988D16

       CALL LC_propagation_checks(L2,C2,nwires-1)
       
!       CALL write_line('FINISHED: LC_checks',0,output_to_screen_flag)
       
       return
       
       end subroutine LC_checks
!
! Name LC_propagation_checks
!     
!
! Description
!     checks the modal propagation characteristics from L and C matrices
!
! Comments:
!      
!
! History
!
!     4/02/15 CJS  Calculate the mode velocities and impedances and write to the cable_info file. 
!			  
!

       SUBROUTINE LC_propagation_checks(L,C,dim)

USE file_information
USE TLM_general

IMPLICIT NONE
       
integer dim
       
real*8 	:: L(dim,dim)
real*8 	:: C(dim,dim)

real*8 	:: CL(dim,dim)

real*8 	:: Tv(dim,dim)
real*8 	:: Ti(dim,dim)

real*8 	:: Tvi(dim,dim)
real*8 	:: Tii(dim,dim)

real*8 	:: ZC(dim,dim)

real*8 	:: Lm(dim,dim)
real*8 	:: Cm(dim,dim)
real*8 	:: Zm(dim)

real*8 	:: D(dim,dim)
real*8 	:: sqrtDI(dim,dim)

real*8 	:: T1(dim,dim),T2(dim,dim),T3(dim,dim)

complex*16	:: gamma(dim)
complex*16	:: Z(dim,dim)
real*8		:: gamma_r(dim)

integer :: row,col,i,j
      
! START       
! CALL write_line('CALLED: LC_propagation_checks',0,output_to_screen_flag)

write(cable_info_file_unit,*)'L(H/m)='
do i=1,dim
  write(cable_info_file_unit,8000)i,(L(i,j),j=1,dim)
end do

write(cable_info_file_unit,*)'C(F/m)='
do i=1,dim
  write(cable_info_file_unit,8000)i,(C(i,j),j=1,dim)
end do

CALL dmatmul(C,dim,dim,L,dim,dim,CL,dim)

! Initially work with CL to obtain CL=[TI]^-1 [ Gamma] [TI]

! calculate the eigenvalues and eigenvectors of CL
CALL deigvectors(CL,dim,GAMMA,Z,dim) 

gamma_r(:)=dble(gamma(:))

write(cable_info_file_unit,*)'Gamma_r'
do row=1,dim
  if (imag(gamma(row)).EQ.0d0) then
    write(cable_info_file_unit,*)row,gamma_r(row)
  else
    write(cable_info_file_unit,*)'****',row,gamma(row),' ****'
  end if
end do

write(cable_info_file_unit,*)'Mode velocities'
do row=1,dim
  write(cable_info_file_unit,*)row,1d0/sqrt(gamma_r(row))
end do

! get the diagonal matrix, D

D(:,:)=0d0
sqrtDI(:,:)=0d0
do row=1,dim
  D(row,row)=gamma_r(row)
  sqrtDI(row,row)=1d0/sqrt(gamma_r(row))
end do

! eigenvector matrix
TI(:,:)=dble(Z(:,:))

CALL dsvd_invert(TI,dim,dim,TII,dim)
write(cable_info_file_unit,*)'Eigenvector matrix:'
write(cable_info_file_unit,*)'TI'
do i=1,dim
  write(cable_info_file_unit,8000)i,(TI(i,j),j=1,dim)
end do

! Calculate the voltage transformation matrices

CALL dtranspose(TII,dim,dim,TV,dim)
CALL dsvd_invert(TV,dim,dim,TVI,dim)

! Calculate the characteristic impedance matrix
CALL dmatmul(sqrtDI,dim,dim,TII,dim,dim,T1,dim)
CALL dmatmul(TI,dim,dim,T1,dim,dim,T2,dim)
CALL dmatmul(L,dim,dim,T2,dim,dim,ZC,dim)

! Calculate the modal inductance matrix
CALL dmatmul(L,dim,dim,TI,dim,dim,T1,dim)
CALL dmatmul(TVI,dim,dim,T1,dim,dim,Lm,dim)

! Calculate the modal capacitance matrix
CALL dmatmul(C,dim,dim,TV,dim,dim,T1,dim)
CALL dmatmul(TII,dim,dim,T1,dim,dim,Cm,dim)

write(cable_info_file_unit,*)'TI_inverse'
do i=1,dim
  write(cable_info_file_unit,8000)i,(TII(i,j),j=1,dim)
end do

write(cable_info_file_unit,*)'TV'
do i=1,dim
  write(cable_info_file_unit,8000)i,(TV(i,j),j=1,dim)
end do

write(cable_info_file_unit,*)'TV_inverse'
do i=1,dim
  write(cable_info_file_unit,8000)i,(TVI(i,j),j=1,dim)
end do

write(cable_info_file_unit,*)'Lm'
do i=1,dim
  write(cable_info_file_unit,8000)i,(Lm(i,j),j=1,dim)
end do

write(cable_info_file_unit,*)'Cm'
do i=1,dim
  write(cable_info_file_unit,8000)i,(Cm(i,j),j=1,dim)
end do

write(cable_info_file_unit,*)'Modal impedances and velocities'
do row=1,dim
  Zm(row)=sqrt(Lm(row,row)/Cm(row,row))
  write(cable_info_file_unit,*)row,Zm(row),1d0/sqrt(Lm(row,row)*Cm(row,row))
end do


! CALL write_line('FINISHED: LC_propagation_checks',0,output_to_screen_flag)
       
  RETURN
     
8000   format (I4,100E14.6)
  
END SUBROUTINE LC_propagation_checks

!
! Name TLM_Cable_checks
!     
!
! Description
!     checks on TLM MTL matrices
!
! Comments:
!      
!
! History
!
!     started 17/03/09 CJS
!
!       
       SUBROUTINE TLM_Cable_checks(Ztl,ZLstub,nwires,flag)

USE file_information
USE TLM_general

IMPLICIT NONE
       
       integer nwires
       real*8 Ztl(nwires,nwires)
       real*8 Ytl(nwires,nwires)
       real*8 Zlstub(nwires,nwires)
       logical flag
       
! local variables            
       real*8 M(nwires,nwires)
       
       complex*16 gamma(nwires)
       integer i,j
       character type

! START
 
     if (flag) then
     
       write(cable_info_file_unit,*)'Ztl:'
       do i=1,nwires
         write(cable_info_file_unit,8000)(Ztl(i,j),j=1,nwires)
       end do
     
       write(cable_info_file_unit,*)'ZLstub:'
       do i=1,nwires
         write(cable_info_file_unit,8000)(ZLstub(i,j),j=1,nwires)
       end do
     end if
     
8000 format (100E14.6)
      
! checks on inductive stub eigenvalues       
       do i=1,nwires
         do j=1,nwires
	   M(i,j)=ZLstub(i,j)
	 end do
       end do
       call deig(M,nwires,GAMMA,nwires)        
       
       if (flag) then
         write(cable_info_file_unit,*)' '
         write(cable_info_file_unit,*)'Eigenvalues of ZLstub'
         do i=1,nwires
           write(cable_info_file_unit,*)i,GAMMA(i)
         end do
         write(cable_info_file_unit,*)' '
       end if
              
       if (flag) then
         write(cable_info_file_unit,*)'Stabilise Ztl'
       end if
       type='T'
       call TLM_stabilise_cable(Ztl,nwires,flag,type)
       
       if (flag) then
         write(cable_info_file_unit,*)'Stabilise ZLstub'
       end if
       type='L'
       call TLM_stabilise_cable(ZLstub,nwires,flag,type)
       
       return
       
       END SUBROUTINE TLM_Cable_checks
!
! Name TLM_stabilise_cable
!     
!
! Description
!     stabilise TLM MTL matrices
!
! Comments:
!      
!
! History
!
!     started 11/01/11 CJS
!
!       
       SUBROUTINE TLM_stabilise_cable(Z,nwires,flag,type)

USE file_information
USE TLM_general

IMPLICIT NONE
       
       integer nwires
       real*8 Z(nwires,nwires)
       logical flag
       character type
       
! local variables       
       
       real*8 M(nwires,nwires)
       real*8 P(nwires,nwires),D(nwires,nwires),PT(nwires,nwires)
       real*8 TM1(nwires,nwires)
       
       real*8 min_positive_eigenvalue
       real*8 max_positive_eigenvalue
       
       integer i,j,row,col
       logical stabilised_flag

! START
       
! copy impedance matrix making it explicitly symmetric      
       do i=1,nwires
         do j=1,nwires
	   M(i,j)=(Z(i,j)+Z(j,i))/2d0
	 end do
       end do
       
       stabilised_flag=.FALSE.
       
       call diagonalise_real_symmetric(M,nwires,P,D,PT,nwires)
       
       if (flag)  then
         write(cable_info_file_unit,*)
         write(cable_info_file_unit,*)'diagonalise'
         write(cable_info_file_unit,*)
       
         write(cable_info_file_unit,*)'M='
         do i=1,nwires
           write(cable_info_file_unit,8000)(M(i,j),j=1,nwires)
         end do
       
         write(cable_info_file_unit,*)'P='
         do i=1,nwires
           write(cable_info_file_unit,8000)(P(i,j),j=1,nwires)
         end do
         write(cable_info_file_unit,*)'GAMMA='
         do i=1,nwires
           write(cable_info_file_unit,8000)D(i,i)
         end do
         write(cable_info_file_unit,*)'PT='
         do i=1,nwires
           write(cable_info_file_unit,8000)(PT(i,j),j=1,nwires)
         end do
8000   format (100E14.6)
       end if
                              
! reform matrix M with positive eigenvalues
      call dmatmul(P,nwires,nwires,D,nwires,nwires,TM1,nwires)
      call dmatmul(TM1,nwires,nwires,PT,nwires,nwires,M,nwires)
       
      if (flag)  then
        write(cable_info_file_unit,*)'Check diagonalisation: M='
        do i=1,nwires
          write(cable_info_file_unit,8000)(M(i,j),j=1,nwires)
        end do
       
        write(cable_info_file_unit,*)'Stabilise eigenvalues:'
	
      end if
      
      min_positive_eigenvalue=1e30
      max_positive_eigenvalue=1d0
      do i=1,nwires
        if (D(i,i).GT.0d0) then
          if (D(i,i).LT.min_positive_eigenvalue) min_positive_eigenvalue=D(i,i)
          if (D(i,i).GT.max_positive_eigenvalue) max_positive_eigenvalue=D(i,i)
        end if
      end do
      do i=1,nwires
        if (D(i,i).LT.0d0) then
          if (type.eq.'L') then
	    D(i,i)=1d-8
            if (flag) write(warning_file_unit,*)'Stabilised eigenvalue in Lstub matrix'
	  else if (type.eq.'T') then
	    D(i,i)=1d-8
            if (flag) write(warning_file_unit,*)'Stabilised eigenvalue in MTL link matrix'
	  else if (type.eq.'C') then
	    D(i,i)=1d8
            if (flag) write(warning_file_unit,*)'Stabilised eigenvalue in Cstub matrix'
	  end if
          stabilised_flag=.TRUE.
        end if
      end do
       
      if (flag)  then
        write(cable_info_file_unit,*)'GAMMA='
        do i=1,nwires
          write(cable_info_file_unit,8000)D(i,i)
        end do
      end if
                              
! reform matrix M with positive eigenvalues
      call dmatmul(P,nwires,nwires,D,nwires,nwires,TM1,nwires)
      call dmatmul(TM1,nwires,nwires,PT,nwires,nwires,M,nwires)
		               
! copy impedance matrix making it explicitly symmetric      
       do i=1,nwires
         do j=1,nwires
	   Z(i,j)=(M(i,j)+M(j,i))/2d0
	 end do
       end do
       
      if (flag)  then
        write(cable_info_file_unit,*)'Reform matrix:'
        do i=1,nwires
          write(cable_info_file_unit,8000)(Z(i,j),j=1,nwires)
        end do
      end if
      
      if (stabilised_flag) then
       
        call diagonalise_real_symmetric(M,nwires,P,D,PT,nwires)
        if (flag)  then
          write(cable_info_file_unit,*)'RECALCULATE GAMMA='
          do i=1,nwires
            write(cable_info_file_unit,8000)D(i,i)
          end do
	
        end if
      end if
              
      return
       
      END SUBROUTINE TLM_stabilise_cable
