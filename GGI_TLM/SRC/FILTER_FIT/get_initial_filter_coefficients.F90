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
! SUBROUTINE get_initial_filter_coefficients
!
! NAME
!     get_initial_filter_coefficients
!
! DESCRIPTION
!
! Fit rational function model to input data using the Wiener Hopf method
! as described in Dawson paper
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 12/12/2012 CJS
!
!
SUBROUTINE get_initial_filter_coefficients

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE filter_types
USE filter_functions
USE constants

IMPLICIT NONE
 
! local variables
! filter definition stuff

  integer :: coeff_dim,freq_dim
  integer :: matdim

  complex*16,allocatable	:: MW(:,:)
  complex*16,allocatable	:: MWI(:,:)
  complex*16,allocatable	:: VBW(:)
  complex*16,allocatable	:: VXW(:)
  
  real*8			:: sigma_average
 
  integer :: i,row,col,function

! START

  CALL write_line('CALLED: get_initial_filter_coefficients',0,ff_output_to_screen)

! Allocate memory
  matdim=2*n_values

  ALLOCATE( MW(matdim,matdim) )
  ALLOCATE( MWI(matdim,matdim) )
  ALLOCATE( VBW(matdim) )
  ALLOCATE( VXW(matdim) )

! allocate filters  
  ALLOCATE( filter_S(1:n_functions) )
  ALLOCATE( filter_S_PR(1:n_functions) )
  ALLOCATE( filter_S_PZ(1:n_functions) )
  ALLOCATE( filter_sigma(1:n_functions) )
  
! loop over number of functions to fit

  write(*,*)'n_values=',n_values
  write(*,*)'n_functions=',n_functions

  do function=1,n_functions
              	      
! fill matrix M
! Include fs at negative frequencies to ensure coefficients are real

    do row=1,n_values

      do i=1,order
    	col=i
    	MW(row,col)=value(function,row)*(s(row)**i)
    	MW(row+n_values,col)=conjg( value(function,row)*(s(row)**i) )
      end do  
      
      do i=0,order
    	col=i+order+1
    	MW(row,col)=-(s(row)**i)
    	MW(row+n_values,col)=-conjg((s(row)**i))
      end do
	 
    end do ! next row

! fill vector V

    do row=1,n_values
      VBW(row)=-value(function,row)
      VBW(row+n_values)=-conjg(value(function,row))
    end do

    freq_dim=2*n_values
    coeff_dim=2*order+1

!    CALL csvd_invert(MW,freq_dim,coeff_dim,MWI,matdim)
!    CALL csvd_invert_LAPACK(MW,freq_dim,coeff_dim,MWI,matdim)

    CALL cinvert_Moore_Penrose(MW,freq_dim,coeff_dim,MWI,matdim) 
    
    CALL cmatvmul(MWI,coeff_dim,freq_dim,VBW,freq_dim,VXW,matdim)

! create filter

    filter_S(function)=allocate_Sfilter(order,order)

    filter_S(function)%wnorm=wnorm
    
    filter_S(function)%b%coeff(0)=1d0
    do i=1,order
      row=i
      filter_S(function)%b%coeff(i)=VXW(row)  
    end do
    do i=0,order
      row=i+order+1
      filter_S(function)%a%coeff(i)=VXW(row)  
    end do

! Calculate the conductivity term if required    
   filter_sigma(function)=0d0   ! assume zero initially
   
   if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 

     sigma_average=0d0
     
     do i=1,n_values
       sigma_average=sigma_average+dble( value(1,i)*s(row) )
     end do
     
     sigma_average=sigma_average/n_values
     
     filter_sigma(function)=sigma_average
     
   end if ! dielectric or magnetic material
   
  end do ! next function
    
! Deallocate memory

  DEALLOCATE( MW )
  DEALLOCATE( MWI )
  DEALLOCATE( VBW )
  DEALLOCATE( VXW )
                   
  RETURN

  CALL write_line('FINISHED: get_initial_filter_coefficients',0,ff_output_to_screen)

END SUBROUTINE get_initial_filter_coefficients
