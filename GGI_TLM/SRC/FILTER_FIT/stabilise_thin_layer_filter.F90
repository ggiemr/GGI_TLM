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
!SUBROUTINE stabilise_thin_layer
!SUBROUTINE correct_Z(Z,Z_correction)
!SUBROUTINE stabilise_thin_layer_OLD
!
! NAME
!     stabilise_thin_layer
!
! DESCRIPTION
!     Given the four impedance filter functions in pole-residue form:
!     1. Ensure that all poles are stable (LHS of s plane)
!     2. Ensure that the impedance matrix is positive definite for all testing frequencies
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/12/2012 CJS
!     5/12/2013 CJS : Improve and simplify the method of applying the correction - the original over-corrected due to an error
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
	
  real*8	:: Z(2,2)
  real*8	:: Z_correction(2,2)

  integer	:: freq_loop

! START

! fill Z matrix with the constant terms from the Pole residue expansion.
! This first stage will ensure passivity of the system at high frequency. 
  Z(1,1)=filter_S_PR(1)%C
  Z(1,2)=filter_S_PR(2)%C
  Z(2,1)=filter_S_PR(3)%C
  Z(2,2)=filter_S_PR(4)%C

! Work out a correction to ensure that Z is positive definite  
  CALL correct_Z(Z,Z_correction)
	 	  
  filter_S_PR(1)%C=filter_S_PR(1)%C+Z_correction(1,1)
  filter_S_PR(2)%C=filter_S_PR(2)%C+Z_correction(1,2)
  filter_S_PR(3)%C=filter_S_PR(3)%C+Z_correction(2,1)
  filter_S_PR(4)%C=filter_S_PR(4)%C+Z_correction(2,2)

! second stage, impose passivity at all frequencies in testing list

! loop over frequencies in testing list	
  
  do freq_loop=1,n_testing_frequencies

! evaluate Zij(s) at this frequency	

    Z(1,1)=dble( evaluate_Sfilter_PR_frequency_response(filter_S_PR(1),testing_frequency(freq_loop)) )
    Z(1,2)=dble( evaluate_Sfilter_PR_frequency_response(filter_S_PR(2),testing_frequency(freq_loop)) )
    Z(2,1)=dble( evaluate_Sfilter_PR_frequency_response(filter_S_PR(3),testing_frequency(freq_loop)) )
    Z(2,2)=dble( evaluate_Sfilter_PR_frequency_response(filter_S_PR(4),testing_frequency(freq_loop)) )

    CALL correct_Z(Z,Z_correction)
	 	  	 	  
    filter_S_PR(1)%C=filter_S_PR(1)%C+Z_correction(1,1)
    filter_S_PR(2)%C=filter_S_PR(2)%C+Z_correction(1,2)
    filter_S_PR(3)%C=filter_S_PR(3)%C+Z_correction(2,1)
    filter_S_PR(4)%C=filter_S_PR(4)%C+Z_correction(2,2)

! next frequency in testing list	   	 	  
  end do

  RETURN

END SUBROUTINE stabilise_thin_layer

! NAME
!     stabilise_thin_layer
!
! DESCRIPTION
!     Given the four impedance filter functions in pole-residue form:
!     1. Ensure that all poles are stable (LHS of s plane)
!     2. Ensure that the impedance matrix is positive definite for all testing frequencies
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/12/2013 CJS
!
!
SUBROUTINE correct_Z(Z,Z_correction)

IMPLICIT NONE
	
  real*8	:: Z(2,2)
  real*8	:: Z_correction(2,2)
  
  real*8	:: Znew(2,2)
  real*8	:: det,trace
  real*8	:: k1,k2,k,delta
 
! local variables

! START

  Z_correction(:,:)=0d0
  
  if ( (Z(1,1).eq.Z(1,2)).AND.(Z(1,1).eq.Z(2,1)).AND.(Z(1,1).eq.Z(2,2)) ) then
    
! We assume that this is a lumped_impedance_model i.e. all impedance parameters are the same
! In this case the determinant will be zero, we only have to ensure that the trace (Z11+Z22)>0 i.e. Z11 and Z22 are both > 0
! We can make this so by replacing any negative values by small positive values

    if (Z(1,1).LT.0d0) then
      Z_correction(:,:)=1.01D0*abs(Z(:,:))
    end if
     
  else
    
    det=Z(1,1)*Z(2,2)-Z(1,2)*Z(2,1)
    trace=Z(1,1)+Z(2,2)

    if ( (det.LT.0d0).OR.(trace.LT.0d0) ) then
! the matrix is not positive definite...

! Start to fix the problem
! 1. make both Z11 and Z22 positive. This will ensure that the trace is positive
      if (Z(1,1).LT.0d0) then
        Z_correction(1,1)=1.01D0*abs(Z(1,1))
      end if
      if (Z(2,2).LT.0d0) then
        Z_correction(2,2)=1.01D0*abs(Z(2,2))
      end if

! the trace should now be positive as both Z11 and Z22 are positive

! check det again     
      Znew=Z+Z_correction
      det=Znew(1,1)*Znew(2,2)-Znew(1,2)*Znew(2,1)
      
      if (det.lt.0d0) then

! multiplication of diagonal values doesn't work too well...      
!        k=1.01d0*Znew(1,2)*Znew(2,1)/(Znew(1,1)*Znew(2,2))	
!	Znew(1,1)=sqrt(k)*Znew(1,1)
!	Znew(2,2)=sqrt(k)*Znew(2,2)

        k1=1.01d0*Znew(1,2)*Znew(2,1)/Znew(2,2)-Znew(1,1)
        k2=1.01d0*Znew(1,2)*Znew(2,1)/Znew(1,1)-Znew(2,2)
	
        if (k1.le.k2) then	
	  Znew(1,1)=Znew(1,1)+k1          
        else
 	  Znew(2,2)=Znew(2,2)+k2     
        end if	
			
      end if
      
      Z_correction=Znew-Z
      
    end if
    
  end if
  
  Znew=Z+Z_correction
  det=Znew(1,1)*Znew(2,2)-Znew(1,2)*Znew(2,1)
  trace=Znew(1,1)+Znew(2,2)

  if ( (det.LT.0d0).OR.(trace.LT.0d0) ) then
    write(*,*)'Error in correct_Z: Could not stabilise the impedance matrix'
    write(*,*)'Orignal Z:'
    write(*,*)Z(1,1),Z(1,2)
    write(*,*)Z(2,1),Z(2,2)
    write(*,*)'Znew:'
    write(*,*)Znew(1,1),Znew(1,2)
    write(*,*)Znew(2,1),Znew(2,2)
    write(*,*)'det=',det
    write(*,*)'trace=',trace
    STOP
  end if
    
  RETURN
  
END SUBROUTINE correct_Z
!
! NAME
!     stabilise_thin_layer_OLD
!
! DESCRIPTION
!     Given the four impedance filter functions in pole-residue form:
!     1. Ensure that all poles are stable (LHS of s plane)
!     2. Ensure that the impedance matrix is positive definite for all testing frequencies
!     
! COMMENTS
!     Trying to improve on this original version which tends to over correct and also has errors...
!
! HISTORY
!
!     started 13/12/2012 CJS
!
!
SUBROUTINE stabilise_thin_layer_OLD

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

  integer	:: i,row,col
  integer	:: freq_loop

  logical	:: lumped_impedance_model

! START

! fill Z matrix with the constant terms from the Pole residue expansion.
! This first stage will ensure passivity of the system at high frequency. 
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
! OLD      correction=correction+1.1D0*abs(GAMMA(row))     
      correction=max( correction,1.1d0*abs(GAMMA(row)) )
    end if
  end do	  

! completed first stage.

  filter_S_PR(1)%C=filter_S_PR(1)%C+correction
  filter_S_PR(4)%C=filter_S_PR(4)%C+correction

! second stage, impose passivity at all frequencies in testing list

! loop over frequencies in testing list	
  
  do freq_loop=1,n_testing_frequencies

! evaluate Zij(s) at this frequency	

    Z(1,1)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(1),testing_frequency(freq_loop))
    Z(1,2)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(2),testing_frequency(freq_loop))
    Z(2,1)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(3),testing_frequency(freq_loop))
    Z(2,2)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(4),testing_frequency(freq_loop))
	  	  
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
	  correction=correction+1.1d0*abs(GAMMA(row))
    	end if
      else
    	if(GAMMA(row).lt.0d0) then
!    	  write(warning_file_unit,*)'Negative eigenvalue found at angular frequency',&
!    		  aimag(stest(i)) ,GAMMA(row),0d0
! OLD	  correction=correction+1.1d0*abs(GAMMA(row))
	  correction=max( correction,1.01d0*abs(GAMMA(row)) )
    	end if
      
      end if
    end do	    

    filter_S_PR(1)%C=filter_S_PR(1)%C+correction
    filter_S_PR(4)%C=filter_S_PR(4)%C+correction

! next frequency in testing list	   	 	  
  end do

  RETURN

END SUBROUTINE stabilise_thin_layer_OLD
