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
  SUBROUTINE allocate_matrix_memory()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  

  integer wire
    
! START
! allocate memory  

! calculate the total number of points and the total number of Fourier terms

  pul_tot_nterms=0
  pul_tot_npoints=0

  do wire=1,pul_nwires
  
    pul_tot_nterms=pul_tot_nterms+1+2*pul_wire_spec(wire)%nterms
    pul_tot_npoints=pul_tot_npoints+pul_wire_spec(wire)%npoints
        
  end do ! next wire
  
  if (pul_include_TLM_return) then
    pul_tot_nterms_dielectric=pul_tot_nterms-1-2*pul_wire_spec(pul_nwires)%nterms
    pul_tot_npoints_dielectric=pul_tot_npoints-pul_wire_spec(pul_nwires)%npoints
  else
    pul_tot_nterms_dielectric=pul_tot_nterms
    pul_tot_npoints_dielectric=pul_tot_npoints
  end if

!  write(*,*)'Total number of conductor potential Fourier terms=',pul_tot_nterms
!  write(*,*)'Total number of dielectric D field Fourier terms =',pul_tot_nterms_dielectric
!  write(*,*)'Total number of conductor potential points=',pul_tot_npoints
!  write(*,*)'Total number of dielectric D field points =',pul_tot_npoints_dielectric
  
  if (pul_tot_npoints+pul_tot_npoints_dielectric.ne.pul_tot_nterms+pul_tot_nterms_dielectric) then
    write(*,*)'Error, total number of potential points should be '
    write(*,*)'equal to the number of Fourier expansion coefficients'
    stop
  end if
  
  pul_LC_matdim=pul_nwires-1
  
  if (pul_include_dielectric) then
    pul_D_matdim=pul_tot_npoints+pul_tot_npoints_dielectric
  else
    pul_D_matdim=pul_tot_npoints
  end if
  
!  write(*,*)'D/B matrix dimension',pul_D_matdim
!  write(*,*)'Cg matrix dimension',pul_nwires
!  write(*,*)'LC matrix dimension',pul_LC_matdim
  
!  write(*,*)'Allocating D/B matrices'
  allocate ( pul_D(1:pul_D_matdim,1:pul_D_matdim) )
  allocate ( pul_B(1:pul_D_matdim,1:pul_D_matdim) )
  
!  write(*,*)'Allocating phi/alpha vectors'
  allocate ( pul_phi(1:pul_D_matdim) )
  allocate ( pul_alpha(1:pul_D_matdim) )
  
!  write(*,*)'Allocating Cg,L,C matrices'
  allocate ( pul_Cg(1:pul_nwires,1:pul_nwires) )
  allocate ( pul_L(1:pul_LC_matdim,1:pul_LC_matdim) )
  allocate ( pul_C(1:pul_LC_matdim,1:pul_LC_matdim) )
  
  RETURN
  
  END SUBROUTINE allocate_matrix_memory
!
! ___________________________________________________________________
!
!   
  SUBROUTINE deallocate_matrix_memory()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  

  integer wire
    
! START
! deallocate matrix memory  

!  write(*,*)'Deallocate D'
  if (allocated ( pul_D ) ) deallocate ( pul_D )
!  write(*,*)'Deallocate B'
  if (allocated ( pul_B ) ) deallocate ( pul_B ) 
!  write(*,*)'Deallocate phi'
  if (allocated ( pul_phi ) ) deallocate ( pul_phi )
!  write(*,*)'Deallocate alpha'
  if (allocated ( pul_alpha ) ) deallocate ( pul_alpha ) 
!  write(*,*)'Deallocate Cg'
  if (allocated ( pul_Cg ) ) deallocate ( pul_Cg )
!  write(*,*)'Deallocate C'
  if (allocated ( pul_C ) ) deallocate ( pul_C )
!  write(*,*)'Deallocate L'
  if (allocated ( pul_L ) ) deallocate ( pul_L )
  
  return
  
  END SUBROUTINE deallocate_matrix_memory
!
! ___________________________________________________________________
!
!   
  SUBROUTINE deallocate_wire_spec_memory()
   
USE pul_data
USE constants

IMPLICIT NONE
  
! local variables  

  integer wire
    
! START
! deallocate memory  
  
!  write(*,*)'Deallocate wires'
  
  do wire=1,pul_nwires
  
    if (allocated( pul_wire_spec(wire)%xp ) ) deallocate ( pul_wire_spec(wire)%xp )
    if (allocated( pul_wire_spec(wire)%yp ) ) deallocate ( pul_wire_spec(wire)%yp )
    if (allocated( pul_wire_spec(wire)%rp ) ) deallocate ( pul_wire_spec(wire)%rp )
    if (allocated( pul_wire_spec(wire)%tp ) ) deallocate ( pul_wire_spec(wire)%tp )
    
    if (allocated( pul_wire_spec(wire)%a ) ) deallocate ( pul_wire_spec(wire)%a )
    if (allocated( pul_wire_spec(wire)%b ) ) deallocate ( pul_wire_spec(wire)%b )
    
    if (pul_include_dielectric) then

      if (allocated( pul_wire_spec(wire)%xdp ) ) deallocate ( pul_wire_spec(wire)%xdp )
      if (allocated( pul_wire_spec(wire)%ydp ) ) deallocate ( pul_wire_spec(wire)%ydp )
      if (allocated( pul_wire_spec(wire)%rdp ) ) deallocate ( pul_wire_spec(wire)%rdp )
      if (allocated( pul_wire_spec(wire)%tdp ) ) deallocate ( pul_wire_spec(wire)%tdp )
      
      if (allocated( pul_wire_spec(wire)%nxdp ) ) deallocate ( pul_wire_spec(wire)%nxdp )
      if (allocated( pul_wire_spec(wire)%nydp ) ) deallocate ( pul_wire_spec(wire)%nydp )
     
      if (allocated( pul_wire_spec(wire)%a2 ) ) deallocate ( pul_wire_spec(wire)%a2 )
      if (allocated( pul_wire_spec(wire)%b2 ) ) deallocate ( pul_wire_spec(wire)%b2 )
      
    end if
    
  end do ! next wire
  
  if (allocated( pul_wire_spec ) ) deallocate ( pul_wire_spec )
  
  return
  
  END SUBROUTINE deallocate_wire_spec_memory
