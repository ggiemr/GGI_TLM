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
!

       SUBROUTINE LC_checks(L,C,length,nwires)

USE file_information
USE TLM_general
USE mesh

IMPLICIT NONE
       
       integer nwires
       real*8 L(nwires,nwires)
       real*8 C(nwires,nwires)
       real*8 length
       
       real*8 LI(nwires,nwires)
       real*8 CI(nwires,nwires)
       real*8 LC(nwires,nwires)
       real*8 LCI(nwires,nwires)
       complex*16 gamma(nwires)
       integer i,j
       logical symmetric_matrix
       
! START       
       CALL write_line('CALLED: LC_checks',0,output_to_screen_flag)

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
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(C(i,j),j=1,nwires)
       end do
       
       write(cable_info_file_unit,*)'Ci='
       write(cable_info_file_unit,*)'Check symmetry Ci'
       call dcheck_symmetry(Ci,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag:',symmetric_matrix
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(Ci(i,j),j=1,nwires)
       end do
       
       write(cable_info_file_unit,*)'L='
       write(cable_info_file_unit,*)'Check symmetry L'
       call dcheck_symmetry(L,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag:',symmetric_matrix
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(L(i,j),j=1,nwires)
       end do
       
       write(cable_info_file_unit,*)'Li='
       write(cable_info_file_unit,*)'Check symmetry Li'
       call dcheck_symmetry(Li,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag:',symmetric_matrix
       do i=1,nwires
         write(cable_info_file_unit,8000)i,(Li(i,j),j=1,nwires)
       end do

8000   format (I4,100E14.6)
       
       write(cable_info_file_unit,*)'calc LC'
       call dmatmul(L,nwires,nwires,C,nwires,nwires,LC,nwires)
       write(cable_info_file_unit,*)'Check symmetry LC'
       call dcheck_symmetry(LC,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag, LC:',symmetric_matrix

       call deig(LC,nwires,GAMMA,nwires)        
       write(cable_info_file_unit,*)' '
       write(cable_info_file_unit,*)'Eigenvalues of LC'
       do i=1,nwires
         write(cable_info_file_unit,*)i,GAMMA(i)
       end do
       write(cable_info_file_unit,*)' '
       write(cable_info_file_unit,*)'Mode velocities'
       do i=1,nwires
         write(cable_info_file_unit,*)i,length/sqrt(GAMMA(i))
       end do
       call dmatmul(L,nwires,nwires,Ci,nwires,nwires,LCi,nwires)
       call dcheck_symmetry(LCi,nwires,nwires,symmetric_matrix,nwires) 
       write(cable_info_file_unit,*)'Symmetry flag, LCi:',symmetric_matrix
       call deig(LCi,nwires,GAMMA,nwires)        
       write(cable_info_file_unit,*)' '
       write(cable_info_file_unit,*)'Eigenvalues of LCi i.e. Mode impedances'
       do i=1,nwires
         write(cable_info_file_unit,*)i,sqrt(GAMMA(i))
       end do
       write(cable_info_file_unit,*)' '
       
       CALL write_line('FINISHED: LC_checks',0,output_to_screen_flag)
       
       return
       
       end subroutine LC_checks

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
