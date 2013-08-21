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
!     stabilise_thin_layer
!
! DESCRIPTION
!     Given the permittivity / permeability filter function in pole-residue form:
!     1. Ensure that all poles are stable (LHS of s plane)
!     2. Ensure that Im{epsr}<=0 for all testing frequencies
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/12/2012 CJS
!
!
SUBROUTINE stabilise_thin_layer

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE FF_file_stuff
USE FF_testing_data
USE filter_functions

IMPLICIT NONE
 
! local variables
	
  complex*16	:: Z(2,2)
  real*8	:: G(2,2),GAMMAM(2,2)
  real*8     	:: GAMMA(2)
  real*8	:: U(2,2),VT(2,2)
  real*8	:: gnorm

  real*8	:: correction
  
  real*8	:: stability_test_small

  integer	:: i,row,col
  integer	:: freq_loop

  logical	:: lumped_impedance_model

! START

  stability_test_small=1D-6

! fill Z matrix
  Z(1,1)=filter_S_PR(1)%C
  Z(1,2)=filter_S_PR(2)%C
  Z(2,1)=filter_S_PR(3)%C
  Z(2,2)=filter_S_PR(4)%C
	 	  
! calculate G as the real part of Z, making G explicitly symmetric
  gnorm=0d0
  do row=1,2
    do col=1,2
      G(row,col)=0.5d0*dble(Z(row,col)+Z(col,row))
      gnorm=gnorm+G(row,col)*G(row,col)
    end do
  end do

  gnorm=sqrt(gnorm)
  				    
  call diagonalise_real_symmetric(G,2,U,GAMMAM,VT,2)
  
  do row=1,2
    GAMMA(row)=GAMMAM(row,row)
  end do
	
! we require all eigenvalues to be positive so check for this
! and make all negative eigenvalues slightly positive (gives
! the optimiser something to work on)	 
  correction=0d0 
  do row=1,2
    if(GAMMA(row).lt.gnorm*stability_test_small) then
!      write(warning_file_unit,*)'Negative eigenvalue found in D matrix'
      correction=correction+1.1D0*abs(GAMMA(row))     
    end if
  end do	  

! completed first stage.

  filter_S_PR(1)%C=filter_S_PR(1)%C+correction
  filter_S_PR(4)%C=filter_S_PR(4)%C+correction

! second stage, impose passivity at all frequencies in testing list

! loop over frequencies in testing list	
  
  do freq_loop=1,n_testing_frequencies

! evaluate Zij(s) at this frequency	

    Z(1,1)=evaluate_Sfilter_frequency_response(filter_S(1),testing_frequency(freq_loop))
    Z(1,2)=evaluate_Sfilter_frequency_response(filter_S(2),testing_frequency(freq_loop))
    Z(2,1)=evaluate_Sfilter_frequency_response(filter_S(3),testing_frequency(freq_loop))
    Z(2,2)=evaluate_Sfilter_frequency_response(filter_S(4),testing_frequency(freq_loop))
	  	  
    if ( (Z(1,1).eq.Z(1,2)).AND.(Z(1,1).eq.Z(2,1)).AND.(Z(1,1).eq.Z(2,2)) ) then
      lumped_impedance_model=.TRUE.
    end if
	 	  
! calculate G as the real part of Z
    gnorm=0.0D0
    do row=1,2
      do col=1,2
    	G(row,col)=dble(Z(row,col))
    	gnorm=gnorm+G(row,col)*G(row,col)
      end do
    end do	    

    gnorm=sqrt(gnorm)

!   write(*,*)'Calling diagonalise_real_symmetric'
    call diagonalise_real_symmetric(G,2,U,GAMMAM,VT,2)

    do row=1,2
      GAMMA(row)=GAMMAM(row,row)
    end do
	  	  
! assemble matrix with only positive eigenvalues 
    correction=0d0
    do row=1,2
      if (.NOT.lumped_impedance_model) then
    	if(GAMMA(row).lt.gnorm*stability_test_small) then
!    	  write(warning_file_unit,*)'Negative eigenvalue found at angular frequency',&
!    		  aimag(stest(i))  ,GAMMA(row),gnorm*stability_test_small
	  correction=correction+1.01d0*abs(GAMMA(row))
    	end if
      else
    	if(GAMMA(row).lt.0d0) then
!    	  write(warning_file_unit,*)'Negative eigenvalue found at angular frequency',&
!    		  aimag(stest(i)) ,GAMMA(row),0d0
	  correction=correction+1.01d0*abs(GAMMA(row))
    	end if
      
      end if
    end do	    

    filter_S_PR(1)%C=filter_S_PR(1)%C+correction
    filter_S_PR(4)%C=filter_S_PR(4)%C+correction

! next frequency in testing list	   	 	  
  end do

  RETURN

END SUBROUTINE stabilise_thin_layer
