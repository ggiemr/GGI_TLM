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
! SUBROUTINE stabilise_FF_input_data
! SUBROUTINE stabilise_FF_input_data_OLD
!
! NAME
!     stabilise_FF_input_data
!
! DESCRIPTION
!     This process removes all data points with unstable (unphysical) data 
!     
! COMMENTS
!     
!    Always make the thin layer reciprocal now
!
! HISTORY
!
!     started 11/12/2012 CJS
!
!
SUBROUTINE stabilise_FF_input_data

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE
 
! local variables

  integer :: freq_loop,function

  complex*16 Z(2,2)
  real*8     G(2,2)
  
  complex*16 GAMMA(2)

  integer row,col
  
  integer			:: new_n_values
  real*8,allocatable		:: new_frequency(:)
  real*8,allocatable		:: new_normalised_frequency(:)
  complex*16,allocatable	:: new_s(:)
  complex*16,allocatable 	:: new_value(:,:)
  
  logical			:: stable_thin_layer

! START

  CALL write_line('CALLED: stabilise_FF_input_data ',0,ff_output_to_screen)

! allocate temporary arrays for stable data
    
  ALLOCATE ( new_frequency(n_values) )
  ALLOCATE ( new_normalised_frequency(n_values) )
  ALLOCATE ( new_s(n_values) )
  ALLOCATE ( new_value(n_functions,n_values) )
  
  new_n_values=0
  
  unstable_input_data=.FALSE.
  
  do freq_loop=1,n_values

    if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 

! ensure that Im{value}.LT.0
      if (Imag(value(1,freq_loop)).GT.0d0) then
      
        unstable_input_data=.TRUE.
	
      else ! stable data at this frequency
      
        new_n_values=new_n_values+1
	new_frequency(new_n_values)=frequency(freq_loop)
	new_normalised_frequency(new_n_values)=normalised_frequency(freq_loop)
	new_s(new_n_values)=s(freq_loop)
	new_value(1,new_n_values)=value(1,freq_loop)
	
      end if
 
    else if (fit_type.eq.thin_layer) then 

! check that Z12=Z21 i.e. that the layer is reciprocal

      if(value(z12,freq_loop).ne.value(z21,freq_loop)) then

! Always make the thin layer reciprocal now
        value(z12,freq_loop)=(value(z12,freq_loop)+value(z21,freq_loop))/2d0
        value(z21,freq_loop)=value(z12,freq_loop)
       
      end if  
     	  
! fill Z matrix 
      Z(1,1)=value(z11,freq_loop)
      Z(1,2)=value(z12,freq_loop)
      Z(2,1)=value(z21,freq_loop)
      Z(2,2)=value(z22,freq_loop)
	 	  
! calculate G as the real part of Z
      
      G(:,:)=dble(Z(:,:))
	  
! calculate eigenvalues of real part of impedance matrix
      call deig(G,2,GAMMA,2)

! we require all eigenvalues to be positive so check for this	  
    
      stable_thin_layer=.TRUE.
      
      do row=1,2
      
        if(dble(GAMMA(row)).lt.small) then
	
          unstable_input_data=.TRUE.
          stable_thin_layer=.FALSE.
	  	
        end if ! GAMMA is small
	
      end do ! next row of GAMMA
      
      if(.NOT.stable_thin_layer) then
	
        unstable_input_data=.TRUE.
	  
      else ! stable data at this frequency
      
        new_n_values=new_n_values+1
        new_frequency(new_n_values)=frequency(freq_loop)
        new_normalised_frequency(new_n_values)=normalised_frequency(freq_loop)
        new_s(new_n_values)=s(freq_loop)
        new_value(1,new_n_values)=value(1,freq_loop)
        new_value(2,new_n_values)=value(2,freq_loop)
        new_value(3,new_n_values)=value(3,freq_loop)
        new_value(4,new_n_values)=value(4,freq_loop)

      end if ! stable thin layer    
     		      
    else if (fit_type.eq.impedance) then 

! ensure that Re{value}.LT.0
      
      if (dble(value(1,freq_loop)).LT.0d0) then
      
        unstable_input_data=.TRUE.
	
      else ! stable data at this frequency
      
        new_n_values=new_n_values+1
	new_frequency(new_n_values)=frequency(freq_loop)
	new_normalised_frequency(new_n_values)=normalised_frequency(freq_loop)
	new_s(new_n_values)=s(freq_loop)
	new_value(1,new_n_values)=value(1,freq_loop)
	
      end if
      
    else if (fit_type.eq.general) then 

! No action to stabilise data 
      
        new_n_values=new_n_values+1
	new_frequency(new_n_values)=frequency(freq_loop)
	new_normalised_frequency(new_n_values)=normalised_frequency(freq_loop)
	new_s(new_n_values)=s(freq_loop)
	new_value(1,new_n_values)=value(1,freq_loop)
   
    end if ! fit_type
    
  end do ! next freq_loop
  
  if (unstable_input_data) then
    CALL write_line('Some of the input data is unstable',0,ff_output_to_screen)
  else
    CALL write_line('The input data is all stable',0,ff_output_to_screen)
  end if 
  
  if (stabilise_input_data_flag) then          
! Copy stable data back to the original array.

    DEALLOCATE ( frequency )
    DEALLOCATE ( normalised_frequency )
    DEALLOCATE ( s )
    DEALLOCATE ( value )
  
    Write(*,*)'Number of input values       =',n_values
    Write(*,*)'Number of stable input values=',new_n_values
  
    n_values=new_n_values  
  
    if (n_values.LT.1) then
      GOTO 9000
    end if
  
    ALLOCATE ( frequency(n_values)  )
    ALLOCATE ( normalised_frequency(n_values)  )
    ALLOCATE ( s(n_values)  )
    ALLOCATE ( value(n_functions,n_values)  )
  
    do freq_loop=1,n_values
    
      frequency(freq_loop)=new_frequency(freq_loop)
      normalised_frequency(freq_loop)=new_normalised_frequency(freq_loop)
      s(freq_loop)=new_s(freq_loop)
      
      do function=1,n_functions
        value(function,freq_loop)=new_value(function,freq_loop)
      end do
      
    end do ! next frequency
    
  end if ! stabilise input data
  
! Deallocate local memory  
  DEALLOCATE ( new_frequency )
  DEALLOCATE ( new_normalised_frequency )
  DEALLOCATE ( new_s )
  DEALLOCATE ( new_value )

  CALL write_line('FINISHED: stabilise_FF_input_data ',0,ff_output_to_screen)

  RETURN
  
9000 CALL write_line('Error in stabilise_FF_input_data',0,.TRUE.)
     CALL write_line('All the input data is unstable',0,.TRUE.)
     STOP

END SUBROUTINE stabilise_FF_input_data
!
! NAME
!     stabilise_FF_input_data_OLD
!
! DESCRIPTION
!     This process makes changes to the input data in order to make it stable...
!     
!     
! COMMENTS
!     
!    Always make the thin layer reciprocal now
!
! HISTORY
!
!     started 11/12/2012 CJS
!
!
SUBROUTINE stabilise_FF_input_data_OLD

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE
 
! local variables

  integer :: freq_loop

  complex*16 Z(2,2)
  real*8     G(2,2)
  
  complex*16 GAMMA(2)

  integer row,col

! START

  CALL write_line('CALLED: stabilise_FF_input_data ',0,ff_output_to_screen)

  unstable_input_data=.FALSE.
  
  do freq_loop=1,n_values

    if ( (fit_type.eq.dielectric_material).OR. (fit_type.eq.magnetic_material)) then 

! ensure that Im{value}.LT.0
      if (Imag(value(1,freq_loop)).GT.0d0) then
      
        unstable_input_data=.TRUE.
	
	if (stabilise_input_data_flag) then
          value(1,freq_loop)=conjg( value(1,freq_loop) )
	end if
	
      end if
 
    else if (fit_type.eq.thin_layer) then 
   
! check that Z12=Z21 i.e. that the layer is reciprocal

      if(value(z12,freq_loop).ne.value(z21,freq_loop)) then

! Always make the thin layer reciprocal now
!        if (stabilise_input_data_flag) then
          value(z12,freq_loop)=(value(z12,freq_loop)+value(z21,freq_loop))/2d0
          value(z21,freq_loop)=value(z12,freq_loop)
!        end if
       
      end if  
     	  
! fill Z matrix 
      Z(1,1)=value(z11,freq_loop)
      Z(1,2)=value(z12,freq_loop)
      Z(2,1)=value(z21,freq_loop)
      Z(2,2)=value(z22,freq_loop)
	 	  
! calculate G as the real part of Z
      
      G(:,:)=dble(Z(:,:))
	  
! calculate eigenvalues of real part of impedance matrix
      call deig(G,2,GAMMA,2)

! we require all eigenvalues to be positive so check for this	  
      do row=1,2
      
        if(dble(GAMMA(row)).lt.small) then
	
          unstable_input_data=.TRUE.
	  
	  if (stabilise_input_data_flag) then
	    GOTO 9000
	  end if
	  
        end if ! GAMMA is small
	
      end do ! next row of GAMMA
     		      
    else if (fit_type.eq.impedance) then 

! ensure that Re{value}.LT.0
      
      if (dble(value(1,freq_loop)).LT.0d0) then
      
        unstable_input_data=.TRUE.
	
	if (stabilise_input_data_flag) then
          value(1,freq_loop)=-conjg( value(1,freq_loop) )
	end if
        
      end if
      
    else if (fit_type.eq.general) then 

! No action
   
    end if ! fit_type
    
  end do ! next freq_loop
  
  if (unstable_input_data) then
    CALL write_line('The input data is unstable',0,ff_output_to_screen)
  else
    CALL write_line('The input data is stable',0,ff_output_to_screen)
  end if 

  CALL write_line('FINISHED: stabilise_FF_input_data ',0,ff_output_to_screen)

  RETURN
  
9000 CALL write_line('Error in stabilise_FF_input_data',0,.TRUE.)
     CALL write_line('Unable to stabilise thin layer data at the moment...',0,.TRUE.)
     STOP

END SUBROUTINE stabilise_FF_input_data_OLD
