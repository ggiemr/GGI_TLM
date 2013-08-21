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
!SUBROUTINE PUL_LC_tlm_interface
!
! NAME
!     SUBROUTINE PUL_LC_tlm_interface
!
! DESCRIPTION
!     Interface to the PUL_LC_CALC process
!     the subroutine takes a list of conductor paramters, creates
!     a bundle cross section geometry then calculates the per unit length L and C 
!     matrices for this cross section with an approriate TLM cell return by calling PUL_LC_CALC
!
!     The TLM return radius is chosen to be different for L and C matrices to try and counteract the
!     propagation velocity error in TLM 
!
!     The PUL parameters are then multiplied by half the TLM cell length to give the per cell segment paramters
!     to create_bundle_LCRG_matrices
!     
! COMMENTS
!     Number of terms in the charge density expansion is fixed at 20 at the moment
!     fictitious_radius_L is set to 1.0/1.08
!     fictitious_radius_C is set to 1.08
!
! HISTORY
!
!     started 20/09/2012 CJS
!
!
SUBROUTINE PUL_LC_tlm_interface(LC_matrix_dimension,conductor_radius,conductor_xc,conductor_yc,	&
                                dielectric_radius,dielectric_permittivity,local_L,local_C,	&
				reference_radius,reference_radius_L,reference_radius_C )

USE TLM_general
USE pul_data
USE constants

IMPLICIT NONE
  
  integer   :: LC_matrix_dimension
  real*8    :: conductor_radius(1:LC_matrix_dimension+1)
  real*8    :: conductor_xc(1:LC_matrix_dimension+1)
  real*8    :: conductor_yc(1:LC_matrix_dimension+1)
  real*8    :: dielectric_radius(1:LC_matrix_dimension+1)
  real*8    :: dielectric_permittivity(1:LC_matrix_dimension+1)
  real*8    :: local_L(1:LC_matrix_dimension,1:LC_matrix_dimension)
  real*8    :: local_C(1:LC_matrix_dimension,1:LC_matrix_dimension)
  
  real*8	:: reference_radius,reference_radius_L,reference_radius_C

! local variables

  real*8	:: fictitious_radius_L
  real*8	:: fictitious_radius_C
  integer	:: wire
  integer	:: tot_n_conductors
  
  integer	:: row,col
  real*8	:: segment_length
  
  real*8	:: cable_radius
  
  real*8,allocatable	:: L_save(:,:)

! START

!  CALL write_line('CALLED: PUL_LC_tlm_interface',0,output_to_screen_flag)

! Set parameters here  

  tot_n_conductors=LC_matrix_dimension+1
  
  if (Cable_LC_Correction_type.EQ.LC_correction_type_geometry_scale) then

! apply the scaling of the fictitious return conductor radius for L and C calculations separately  
    fictitious_radius_L=Inductance_equivalent_radius_factor
    fictitious_radius_C=Capacitance_equivalent_radius_factor
    
  else

! Set the fictitious return conductor radius for L and C calculations to be 1.0
    fictitious_radius_L=1d0
    fictitious_radius_C=1d0
  
  end if
  
  cable_radius=dielectric_radius(LC_matrix_dimension+1)
  
  reference_radius=(dl*TLM_cell_equivalent_radius_factor/2d0)
  
  if (cable_radius.GT.reference_radius*Max_cable_bundle_diameter_factor) then 
    reference_radius=cable_radius/Max_cable_bundle_diameter_factor
  end if  
  
  reference_radius_L=reference_radius*fictitious_radius_L
  reference_radius_C=reference_radius*fictitious_radius_C
    
  pul_tot_nterms=number_of_Fourier_terms_in_PUL_LC_calc  ! number of Fourier series terms for charge density expansion
  
! put the bundle geometry into the format of pul_LC_calc... 
     
  pul_dl=dl
  pul_include_TLM_return=.TRUE.      
  pul_nwires_in=LC_matrix_dimension
  pul_nwires=pul_nwires_in+1  ! include TLM return conductor
  
  pul_return_conductor=pul_nwires
  
  pul_op_flag=.FALSE.
      
  allocate ( pul_wire_spec(1:pul_nwires) )

  do wire=1,pul_nwires_in
  
    pul_wire_spec(wire)%nterms  =pul_tot_nterms
    pul_wire_spec(wire)%npoints=1+2*pul_wire_spec(wire)%nterms
    pul_wire_spec(wire)%xc  =conductor_xc(wire)
    pul_wire_spec(wire)%yc  =conductor_yc(wire)
    pul_wire_spec(wire)%rw  =conductor_radius(wire)
    pul_wire_spec(wire)%ri  =dielectric_radius(wire)
    pul_wire_spec(wire)%epsr=dielectric_permittivity(wire)
        	
  end do ! next wire

! include TLM return  
  wire=pul_nwires
  pul_wire_spec(wire)%xc=0d0
  pul_wire_spec(wire)%yc=0d0
      
  pul_wire_spec(wire)%rw=reference_radius_L ! fictitious radius for inductance calculation      
  pul_wire_spec(wire)%ri=pul_wire_spec(wire)%rw
  pul_wire_spec(wire)%epsr=1.0
  pul_wire_spec(wire)%nterms=pul_tot_nterms
  pul_wire_spec(wire)%npoints=1+2*pul_wire_spec(wire)%nterms
  
  conductor_radius(wire)=pul_wire_spec(wire)%rw
  dielectric_radius(wire)=pul_wire_spec(wire)%ri
  dielectric_permittivity(wire)=pul_wire_spec(wire)%epsr

! This process already done in create_bundle_LCRG_matrices      
!! work out the geometric cross section of the bundle  
!  
!  CALL create_bundle_cross_section_geometry(	LC_matrix_dimension,		&
!      						conductor_radius,dielectric_radius,	&
!      						conductor_xc,conductor_yc)
 
! copy conductor centre coordinates to pul_wire_spec array      

  do wire=1,pul_nwires
      
    pul_wire_spec(wire)%xc=conductor_xc(wire)
    pul_wire_spec(wire)%yc=conductor_yc(wire)
    pul_wire_spec(wire)%rw=conductor_radius(wire)
    pul_wire_spec(wire)%ri=dielectric_radius(wire)
    
    write(*,*)'wire',wire,' xc=',conductor_xc(wire),' yc=',conductor_yc(wire)
        	
  end do ! next wire      

! solve for per-unit length parameters

! solve for capacitance matrix without dielectric first - this then leads to the inductance matrix
  pul_include_dielectric=.FALSE. 
  pul_include_TLM_return=.TRUE.
  pul_op_flag=.FALSE.
                  
  call calc_conductor_points()
  
  call calc_dielectric_points()
  
  call allocate_matrix_memory()  
  
  call calc_D_matrix()
    
! Invert D matrix 
  call dsvd_invert(pul_D,pul_D_matdim,pul_D_matdim,pul_B,pul_D_matdim) 
  
! calculate generalised capacitance matrix  
  call calc_Cg_matrix() 
    
! calculate capacitance matrix  
  call calc_C_matrix()
      
! Calculate the Inductance matrix. L=mu. eps C^-1
  
  call dsvd_invert(pul_C,pul_LC_matdim,pul_LC_matdim,pul_L,pul_LC_matdim)
    
  pul_L(:,:)=pul_L(:,:)*mu0*eps0  
  
!  write(*,*)'L(uH/m):'
!  call write_matrix(pul_L,pul_LC_matdim,1e6,0)
      
  allocate( L_save(1:pul_LC_matdim,1:pul_LC_matdim) )
  L_save(:,:)=pul_L(:,:)
       
! We must recalculate C including the dielectric

!  write(*,*)'Deallocate matrix memory'
  call deallocate_matrix_memory()
        
  pul_include_dielectric=.TRUE.
      
  wire=pul_nwires
      
  pul_wire_spec(wire)%rw=reference_radius_C ! fictitious radius for capacitance calculation            
  pul_wire_spec(wire)%ri=pul_wire_spec(wire)%rw
  pul_wire_spec(wire)%epsr=1.0
  pul_wire_spec(wire)%nterms=pul_tot_nterms
  pul_wire_spec(wire)%npoints=1+2*pul_wire_spec(wire)%nterms
      	  
  call calc_conductor_points()
  
  call calc_dielectric_points()
	
  call allocate_matrix_memory()
  
! fill matrices  
  call calc_D_matrix_dielectric()
  
! Invert D matrix  
  call dsvd_invert(pul_D,pul_D_matdim,pul_D_matdim,pul_B,pul_D_matdim)

! calculate generalised capacitance matrix  
  call calc_Cg_matrix()
    
! calculate capacitance matrix  
  call calc_C_matrix()
      
  segment_length=1d0  ! return per-unit length L and C from here
      
  do row=1,LC_matrix_dimension 
    do col=1,LC_matrix_dimension 
    
      local_L(row,col)=L_save(row,col)*segment_length
      local_C(row,col)=pul_C(row,col) *segment_length

    end do ! next column
  end do ! next row

  deallocate( L_save )

  call deallocate_matrix_memory()
  call deallocate_wire_spec_memory()

!  CALL write_line('FINISHED: PUL_LC_tlm_interface',0,output_to_screen_flag)
    
  RETURN
  
  
END SUBROUTINE PUL_LC_tlm_interface
