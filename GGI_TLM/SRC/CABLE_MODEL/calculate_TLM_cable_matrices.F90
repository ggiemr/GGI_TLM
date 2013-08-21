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
!     calculate_TLM_cable_matrices
!
! DESCRIPTION
!     
! Set cable link and stub matrices for the TLM cable update     
!
! COMMENTS
!     
! Stub calculation: - we need to decide what we need to save on the segment and what is just
! a local array required for the calculation...
!
!
! HISTORY
!
!     started 21/09/2012 CJS
!
!
SUBROUTINE calculate_TLM_cable_matrices(L,C,Ztl,Zlstub,n_conductors)

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE Cables
USE file_information

IMPLICIT NONE

  integer   :: n_conductors
  real*8    :: L(n_conductors,n_conductors)
  real*8    :: C(n_conductors,n_conductors)
  real*8    :: Ztl(n_conductors,n_conductors)
  real*8    :: ZLstub(n_conductors,n_conductors)

! local variables

  real*8    :: Ci(n_conductors,n_conductors)
    
  real*8    :: dt2

! START

!  CALL write_line('CALLED: calculate_TLM_cable_matrices',0,output_to_screen_flag)

  dt2=dt/2d0

! Stub calculation 

  call dsvd_invert(C,n_conductors,n_conductors,Ci,n_conductors)
		    
  Ztl(:,:)=(Ci(:,:)*dt2)  ! include all of the capacitance - no capacitive stubs needed here
		     
  ZLstub(:,:)=2.0*L(:,:)/dt-Ztl(:,:)
  
!  CALL write_line('FINISHED: calculate_TLM_cable_matrices',0,output_to_screen_flag)

  RETURN

END SUBROUTINE calculate_TLM_cable_matrices
