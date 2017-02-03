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
! SUBROUTINE test_stability
! SUBROUTINE test_stability_PZ
! SUBROUTINE test_stability_PR
!
! NAME
!     test_stability
!
! DESCRIPTION
!     Evaluate the filter function at all testing frequencies and test for stability of the model
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE test_stability(write_error_to_file)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_testing_data
USE FF_file_stuff
USE filter_functions
USE constants

IMPLICIT NONE

  logical	:: write_error_to_file
 
! local variables

  integer 	:: function_loop,i
  
  integer 	:: freq_loop

  complex*16 	:: Z(2,2)
  real*8 	:: G(2,2)
  
  complex*16 	:: GAMMA(2)
  
  complex*16	:: Z2
  complex*16	:: epsr_mur

  integer row,col

! START

  CALL write_line('CALLED: test_stability ',0,ff_output_to_screen)

! Test whether the poles of the filter(s) are on the LHS of the S plane                  
  do function_loop=1,n_functions
    filter_S_PR(function_loop)=Convert_filter_S_to_S_PR( filter_S(function_loop) )
    do i=1,filter_S_PR(function_loop)%order
      if (dble(filter_S_PR(function_loop)%poles(i)).gt.(-stability_test_small)) then ! was .GT.0d0
        stable_filter=.FALSE.
      end if
    end do
  end do 
  
! Test whether the conductivity is negative
  if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 
    if (filter_sigma(1).LT.0d0) stable_filter=.FALSE.
  end if
     
  do freq_loop=1,n_testing_frequencies

    if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 
    
      epsr_mur=evaluate_Sfilter_frequency_response(filter_S(1),testing_frequency(freq_loop))
      
! ensure that Im{epsr_mur}.LT.0
      if (Imag(epsr_mur).GT.0d0) then
      
        stable_filter=.FALSE.
	
      end if

! test high frequency value GT 1      
      if ( (freq_loop.EQ.n_testing_frequencies).AND.(dble(epsr_mur).LT.1d0) ) then
      
        stable_filter=.FALSE.
	
      end if
      
    else if (fit_type.eq.thin_layer) then 
        	  
! fill Z matrix 
      Z(1,1)=evaluate_Sfilter_frequency_response(filter_S(1),testing_frequency(freq_loop))
      Z(1,2)=evaluate_Sfilter_frequency_response(filter_S(2),testing_frequency(freq_loop))
      Z(2,1)=evaluate_Sfilter_frequency_response(filter_S(3),testing_frequency(freq_loop))
      Z(2,2)=evaluate_Sfilter_frequency_response(filter_S(4),testing_frequency(freq_loop))
	 	  
! calculate G as the real part of Z
      
      G(:,:)=dble(Z(:,:))
	  
! calculate eigenvalues of real part of impedance matrix
      call deig(G,2,GAMMA,2)

! we require all eigenvalues to be positive so check for this	  
      do row=1,2
      
        if(dble(GAMMA(row)).lt.stability_test_small) then
	        
          stable_filter=.FALSE.
	  
        end if ! GAMMA is small
	
      end do ! next row of GAMMA
     		      
    else if (fit_type.eq.impedance) then 
    
      Z2=evaluate_Sfilter_frequency_response(filter_S(1),testing_frequency(freq_loop))

! ensure that Re{Z2}.LT.0
      
      if (dble(Z2).LT.0d0) then
           
        stable_filter=.FALSE.
        
      end if
   
     		      
    else if (fit_type.eq.general) then 
    
! No action
   
    end if ! fit_type
    
  end do ! next freq_loop
  
  if (.NOT.stable_filter) then
    CALL write_line('The filter is UNSTABLE',0,ff_output_to_screen)
    if (write_error_to_file) then
      open(UNIT=local_file_unit,FILE='Filter_fit_error.dat')
      write(local_file_unit,8000)'Final Sfilter, order',order,' Mean square error=',Mean_square_error,': UNSTABLE '
8000 format(A20,I5,A19,F12.4,A11)
      close(UNIT=local_file_unit)
    end if
  else
    CALL write_line('The filter is STABLE',0,ff_output_to_screen)
    if (write_error_to_file) then
      open(UNIT=local_file_unit,FILE='Filter_fit_error.dat')
      write(local_file_unit,8000)'Final Sfilter, order',order,' Mean square error=',Mean_square_error,':   STABLE'
      close(UNIT=local_file_unit)
    end if
  end if 

  CALL write_line('FINISHED: test_stability ',0,ff_output_to_screen)

  RETURN

END SUBROUTINE test_stability
!
! NAME
!     test_stability_PZ
!
! DESCRIPTION
!     Evaluate the filter function at all testing frequencies and test for stability of the model
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE test_stability_PZ

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_testing_data
USE FF_file_stuff
USE filter_functions
USE constants

IMPLICIT NONE
 
! local variables

  integer 	:: function_loop,i
  
  integer 	:: freq_loop

  complex*16 	:: Z(2,2)
  real*8 	::     G(2,2)
  
  complex*16 	:: GAMMA(2)
  
  complex*16	:: Z2
  complex*16	:: epsr_mur

  integer row,col

! START

!  CALL write_line('CALLED: test_stability_PZ ',0,ff_output_to_screen)

  stable_filter=.TRUE.

! Check the stability of the poles i.e. are they on the LHS of the S plane?                  
  do function_loop=1,n_functions
    do i=1,filter_S_PZ(function_loop)%order
      if (dble(filter_S_PZ(function_loop)%poles(i)).gt.(-stability_test_small)) then ! was .GT.0d0
        stable_filter=.FALSE.
      end if
    end do
  end do 
  
! Test whether the conductivity is negative
  if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 
    if (filter_sigma(1).LT.0d0) stable_filter=.FALSE.
  end if
   
  do freq_loop=1,n_testing_frequencies

    if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 
    
      epsr_mur=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(1),testing_frequency(freq_loop))
      
! ensure that Im{epsr_mur}.LT.0
      if (Imag(epsr_mur).GT.0d0) then
      
        stable_filter=.FALSE.
	
      end if

! test high frequency value GT 1      
      if ( (freq_loop.EQ.n_testing_frequencies).AND.(dble(epsr_mur).LT.1d0) ) then
      
        stable_filter=.FALSE.
	
      end if
 
    else if (fit_type.eq.thin_layer) then 
        	  
! fill Z matrix 
      Z(1,1)=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(1),testing_frequency(freq_loop))
      Z(1,2)=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(2),testing_frequency(freq_loop))
      Z(2,1)=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(3),testing_frequency(freq_loop))
      Z(2,2)=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(4),testing_frequency(freq_loop))
	 	  
! calculate G as the real part of Z
      
      G(:,:)=dble(Z(:,:))
	  
! calculate eigenvalues of real part of impedance matrix
      call deig(G,2,GAMMA,2)

! we require all eigenvalues to be positive so check for this	  
      do row=1,2
      
        if(dble(GAMMA(row)).lt.stability_test_small) then
	        
          stable_filter=.FALSE.
	  
        end if ! GAMMA is small
	
      end do ! next row of GAMMA
     		      
    else if (fit_type.eq.impedance) then 
    
      Z2=evaluate_Sfilter_PZ_frequency_response(filter_S_PZ(1),testing_frequency(freq_loop))

! ensure that Re{Z2}.LT.0
      
      if (dble(Z2).LT.0d0) then
           
        stable_filter=.FALSE.
        
      end if
     		      
    else if (fit_type.eq.general) then 
   
! No action
   
    end if ! fit_type
    
  end do ! next freq_loop
  
!  if (.NOT.stable_filter) then
!    CALL write_line('The filter is UNSTABLE',0,ff_output_to_screen)
!  else
!    CALL write_line('The filter is STABLE',0,ff_output_to_screen)
!  end if 

!  CALL write_line('FINISHED: test_stability_PZ ',0,ff_output_to_screen)

  RETURN

END SUBROUTINE test_stability_PZ
!
! NAME
!     test_stability_PR
!
! DESCRIPTION
!     Evaluate the filter function at all testing frequencies and test for stability of the model
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/12/2013 CJS
!
!
SUBROUTINE test_stability_PR

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_testing_data
USE FF_file_stuff
USE filter_functions
USE constants

IMPLICIT NONE
 
! local variables

  integer 	:: function_loop,i
  
  integer 	:: freq_loop

  complex*16 	:: Z(2,2)
  real*8 	::     G(2,2)
  
  complex*16 	:: GAMMA(2)
  
  complex*16	:: Z2
  complex*16	:: epsr_mur

  integer row,col

! START

!  CALL write_line('CALLED: test_stability_PR ',0,ff_output_to_screen)

  stable_filter=.TRUE.

! Check the stability of the poles i.e. are they on the LHS of the S plane?                  
  do function_loop=1,n_functions
    do i=1,filter_S_PR(function_loop)%order
      if (dble(filter_S_PR(function_loop)%poles(i)).gt.(-stability_test_small)) then ! was .GT.0d0
        stable_filter=.FALSE.
      end if
    end do
  end do 
  
! Test whether the conductivity is negative
  if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 
    if (filter_sigma(1).LT.0d0) stable_filter=.FALSE.
  end if
  
  do freq_loop=1,n_testing_frequencies

    if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 
    
      epsr_mur=evaluate_Sfilter_PR_frequency_response(filter_S_PR(1),testing_frequency(freq_loop))
      
! ensure that Im{epsr_mur}.LT.0
      if (Imag(epsr_mur).GT.0d0) then
      
        stable_filter=.FALSE.
	
      end if

! test high frequency value GT 1      
      if ( (freq_loop.EQ.n_testing_frequencies).AND.(dble(epsr_mur).LT.1d0) ) then
      
        stable_filter=.FALSE.
	
      end if
 
    else if (fit_type.eq.thin_layer) then 
        	  
! fill Z matrix 
      Z(1,1)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(1),testing_frequency(freq_loop))
      Z(1,2)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(2),testing_frequency(freq_loop))
      Z(2,1)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(3),testing_frequency(freq_loop))
      Z(2,2)=evaluate_Sfilter_PR_frequency_response(filter_S_PR(4),testing_frequency(freq_loop))
	 	  
! calculate G as the real part of Z
      
      G(:,:)=dble(Z(:,:))
	  
! calculate eigenvalues of real part of impedance matrix
      call deig(G,2,GAMMA,2)

! we require all eigenvalues to be positive so check for this	  
      do row=1,2
      
        if(dble(GAMMA(row)).lt.stability_test_small) then
	        
          stable_filter=.FALSE.
	  
        end if ! GAMMA is small
	
      end do ! next row of GAMMA
     		      
    else if (fit_type.eq.impedance) then 
    
      Z2=evaluate_Sfilter_PR_frequency_response(filter_S_PR(1),testing_frequency(freq_loop))

! ensure that Re{Z2}.LT.0
      
      if (dble(Z2).LT.0d0) then
           
        stable_filter=.FALSE.
        
      end if
     		      
    else if (fit_type.eq.general) then 
   
! No action
   
    end if ! fit_type
    
  end do ! next freq_loop
  
!  if (.NOT.stable_filter) then
!    CALL write_line('The filter is UNSTABLE',0,ff_output_to_screen)
!  else
!    CALL write_line('The filter is STABLE',0,ff_output_to_screen)
!  end if 

!  CALL write_line('FINISHED: test_stability_PR ',0,ff_output_to_screen)

  RETURN

END SUBROUTINE test_stability_PR
