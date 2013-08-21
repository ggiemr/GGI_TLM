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
! SUBROUTINE dielectric_magnetic_S_PZ_filter_to_opt_param_list
! SUBROUTINE opt_param_list_to_dielectric_magnetic_S_PZ_filter(opt_dim,params)
! SUBROUTINE thin_layer_S_PZ_filter_to_opt_param_list
! SUBROUTINE opt_param_list_to_thin_layer_S_PZ_filter(opt_dim,params)
! SUBROUTINE impedance_S_PZ_filter_to_opt_param_list
! SUBROUTINE opt_param_list_to_impedance_S_PZ_filter(opt_dim,params)
!
! NAME
!     dielectric_magnetic_S_PZ_filter_to_opt_param_list
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE dielectric_magnetic_S_PZ_filter_to_opt_param_list(opt_dim,params)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE

  integer		:: opt_dim

  real*8		:: params(1:opt_dim) 
 
! local variables

  integer		:: i
  integer		:: param_count
  integer		:: pole_count
  integer		:: zero_count

! START 
 
! check...
  if (opt_dim.NE.(order*2+2)) then
    write(*,*)'Error in dielectric_magnetic_S_PZ_filter_to_opt_param_list'
    write(*,*)'opt_dim discrepancy',opt_dim,order*2+2
    STOP
  end if
  
  param_count=0

! Gain term
  param_count=param_count+1
  params(param_count)=filter_S_PZ(1)%G
  
! Poles  

  pole_count=0
  
! real poles first

  do i=1,filter_S_PZ(1)%n_real_poles
    param_count=param_count+1
    pole_count=pole_count+1
    params(param_count)=dble(filter_S_PZ(1)%poles(pole_count))
  end do
  
! complex poles second

  do i=1,filter_S_PZ(1)%n_complex_pole_pairs
    param_count=param_count+1
    pole_count=pole_count+1
    params(param_count)=dble(filter_S_PZ(1)%poles(pole_count))
 ! conjugate pole    
    param_count=param_count+1
    pole_count=pole_count+1
    params(param_count)=imag(filter_S_PZ(1)%poles(pole_count))
  end do
  
! Zeros  

  zero_count=0
  
! real zeros first

  do i=1,filter_S_PZ(1)%n_real_zeros
    param_count=param_count+1
    zero_count=zero_count+1
    params(param_count)=dble(filter_S_PZ(1)%zeros(zero_count))
  end do
  
! complex zeros second

  do i=1,filter_S_PZ(1)%n_complex_zero_pairs
    param_count=param_count+1
    zero_count=zero_count+1
    params(param_count)=dble(filter_S_PZ(1)%zeros(zero_count))
 ! conjugate zero    
    param_count=param_count+1
    zero_count=zero_count+1
    params(param_count)=imag(filter_S_PZ(1)%zeros(zero_count))
  end do
  
! Conductivity
  param_count=param_count+1
  params(param_count)=filter_sigma(1)
  
! check we have set the expected number of parameters
  if (param_count.NE.opt_dim) then
    write(*,*)'Error in dielectric_magnetic_S_PZ_filter_to_opt_param_list'
    write(*,*)'opt_dim discrepancy',opt_dim,param_count
    STOP 
  end if
       
END SUBROUTINE dielectric_magnetic_S_PZ_filter_to_opt_param_list
!
! NAME
!     opt_param_list_to_dielectric_magnetic_S_PZ_filter
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE opt_param_list_to_dielectric_magnetic_S_PZ_filter(opt_dim,params)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE

  integer		:: opt_dim

  real*8		:: params(1:opt_dim) 
 
! local variables

  integer		:: i
  integer		:: param_count
  integer		:: pole_count
  integer		:: zero_count

! START 

! check...
  if (opt_dim.NE.(order*2+2)) then
    write(*,*)'Error in opt_param_list_to_dielectric_magnetic_S_PZ_filter'
    write(*,*)'opt_dim discrepancy',opt_dim,order*2+2
    STOP
  end if
  
  param_count=0

! Gain term
  param_count=param_count+1
  filter_S_PZ(1)%G=params(param_count)
  
! Poles  

  pole_count=0
  
! real poles first

  do i=1,filter_S_PZ(1)%n_real_poles
    param_count=param_count+1
    pole_count=pole_count+1
    filter_S_PZ(1)%poles(pole_count)=params(param_count)
  end do
  
! complex poles second

  do i=1,filter_S_PZ(1)%n_complex_pole_pairs
    param_count=param_count+1
    pole_count=pole_count+1
    filter_S_PZ(1)%poles(pole_count)  =DCMPLX(params(param_count), params(param_count+1))
    filter_S_PZ(1)%poles(pole_count+1)=DCMPLX(params(param_count),-params(param_count+1))
    param_count=param_count+1
    pole_count=pole_count+1
  end do
  
! Zeros  

  zero_count=0
  
! real zeros first

  do i=1,filter_S_PZ(1)%n_real_zeros
    param_count=param_count+1
    zero_count=zero_count+1
    filter_S_PZ(1)%zeros(zero_count)=params(param_count)
  end do
  
! complex zeros second

  do i=1,filter_S_PZ(1)%n_complex_zero_pairs
    param_count=param_count+1
    zero_count=zero_count+1
    filter_S_PZ(1)%zeros(zero_count)  =DCMPLX(params(param_count), params(param_count+1))
    filter_S_PZ(1)%zeros(zero_count+1)=DCMPLX(params(param_count),-params(param_count+1))
    param_count=param_count+1
    zero_count=zero_count+1
  end do
  
! Conductivity
  param_count=param_count+1
  filter_sigma(1)=params(param_count)
  
! check we have set the expected number of parameters
  if (param_count.NE.opt_dim) then
    write(*,*)'Error in opt_param_list_to_dielectric_magnetic_S_PZ_filter'
    write(*,*)'opt_dim discrepancy',opt_dim,param_count
    STOP 
  end if
       
END SUBROUTINE opt_param_list_to_dielectric_magnetic_S_PZ_filter
!
! NAME
!     thin_layer_S_PZ_filter_to_opt_param_list
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE thin_layer_S_PZ_filter_to_opt_param_list(opt_dim,params)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE

  integer		:: opt_dim

  real*8		:: params(1:opt_dim) 
 
! local variables

  integer		:: loop
  integer		:: loop_to_filter_function(1:3)
  integer		:: filter_function
  integer		:: i
  integer		:: param_count
  integer		:: pole_count
  integer		:: zero_count

! START 

! check...
  if (opt_dim.NE.3*(order*2+1)) then
    write(*,*)'Error in thin_layer_S_PZ_filter_to_opt_param_list'
    write(*,*)'opt_dim discrepancy',opt_dim,3*(order*2+1)
    STOP
  end if
  
  param_count=0
  
  loop_to_filter_function(1)=1
  loop_to_filter_function(2)=2
  loop_to_filter_function(3)=4
  
  do loop=1,3
  
    filter_function=loop_to_filter_function(loop)
    
! Gain term
    param_count=param_count+1
    params(param_count)=filter_S_PZ(filter_function)%G
  
! Poles  

    pole_count=0
  
! real poles first

    do i=1,filter_S_PZ(filter_function)%n_real_poles
      param_count=param_count+1
      pole_count=pole_count+1
      params(param_count)=dble(filter_S_PZ(filter_function)%poles(pole_count))
    end do
  
! complex poles second

    do i=1,filter_S_PZ(filter_function)%n_complex_pole_pairs
      param_count=param_count+1
      pole_count=pole_count+1
      params(param_count)=dble(filter_S_PZ(filter_function)%poles(pole_count))
 ! conjugate pole    
      param_count=param_count+1
      pole_count=pole_count+1
      params(param_count)=imag(filter_S_PZ(filter_function)%poles(pole_count))
    end do
  
! Zeros  

    zero_count=0
  
! real zeros first

    do i=1,filter_S_PZ(filter_function)%n_real_zeros
      param_count=param_count+1
      zero_count=zero_count+1
      params(param_count)=dble(filter_S_PZ(filter_function)%zeros(zero_count))
    end do
  
! complex zeros second

    do i=1,filter_S_PZ(filter_function)%n_complex_zero_pairs
      param_count=param_count+1
      zero_count=zero_count+1
      params(param_count)=dble(filter_S_PZ(filter_function)%zeros(zero_count))
 ! conjugate zero    
      param_count=param_count+1
      zero_count=zero_count+1
      params(param_count)=imag(filter_S_PZ(filter_function)%zeros(zero_count))
    end do
     
  end do ! next filter function   
    
! check we have set the expected number of parameters
  if (param_count.NE.opt_dim) then
    write(*,*)'Error in thin_layer_S_PZ_filter_to_opt_param_list'
    write(*,*)'opt_dim discrepancy',opt_dim,param_count
    STOP 
  end if
       
END SUBROUTINE thin_layer_S_PZ_filter_to_opt_param_list
!
! NAME
!     opt_param_list_to_thin_layer_S_PZ_filter
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE opt_param_list_to_thin_layer_S_PZ_filter(opt_dim,params)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE

  integer		:: opt_dim

  real*8		:: params(1:opt_dim) 
 
! local variables
  integer		:: loop
  integer		:: loop_to_filter_function(1:3)
  integer		:: filter_function

  integer		:: i
  integer		:: param_count
  integer		:: pole_count
  integer		:: zero_count

! START 

! check...
  if (opt_dim.NE.3*(order*2+1)) then
    write(*,*)'Error in opt_param_list_to_thin_layer_S_PZ_filter'
    write(*,*)'opt_dim discrepancy',opt_dim,3*(order*2+1)
    STOP
  end if
  
  param_count=0
  
  loop_to_filter_function(1)=1
  loop_to_filter_function(2)=2
  loop_to_filter_function(3)=4
  
  do loop=1,3
  
    filter_function=loop_to_filter_function(loop)

! Gain term
    param_count=param_count+1
    filter_S_PZ(filter_function)%G=params(param_count)
  
! Poles  

    pole_count=0
  
! real poles first

    do i=1,filter_S_PZ(filter_function)%n_real_poles
      param_count=param_count+1
      pole_count=pole_count+1
      filter_S_PZ(filter_function)%poles(pole_count)=params(param_count)
    end do
  
! complex poles second

    do i=1,filter_S_PZ(filter_function)%n_complex_pole_pairs
      param_count=param_count+1
      pole_count=pole_count+1
      filter_S_PZ(filter_function)%poles(pole_count)  =DCMPLX(params(param_count), params(param_count+1))
      filter_S_PZ(filter_function)%poles(pole_count+1)=DCMPLX(params(param_count),-params(param_count+1))
      param_count=param_count+1
      pole_count=pole_count+1
    end do
  
! Zeros  

    zero_count=0
  
! real zeros first

    do i=1,filter_S_PZ(filter_function)%n_real_zeros
      param_count=param_count+1
      zero_count=zero_count+1
      filter_S_PZ(filter_function)%zeros(zero_count)=params(param_count)
    end do
  
! complex zeros second

    do i=1,filter_S_PZ(filter_function)%n_complex_zero_pairs
      param_count=param_count+1
      zero_count=zero_count+1
      filter_S_PZ(filter_function)%zeros(zero_count)  =DCMPLX(params(param_count), params(param_count+1))
      filter_S_PZ(filter_function)%zeros(zero_count+1)=DCMPLX(params(param_count),-params(param_count+1))
      param_count=param_count+1
      zero_count=zero_count+1
    end do
  
  end do ! next filter function
  
! Reciprocal impedance matrix do Z21=Z12
  filter_S_PZ(3)=filter_S_PZ(2)
  
! check we have set the expected number of parameters
  if (param_count.NE.opt_dim) then
    write(*,*)'Error in opt_param_list_to_thin_layer_S_PZ_filter'
    write(*,*)'opt_dim discrepancy',opt_dim,param_count
    STOP 
  end if
       
END SUBROUTINE opt_param_list_to_thin_layer_S_PZ_filter
!
! NAME
!     impedance_S_PZ_filter_to_opt_param_list
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE impedance_S_PZ_filter_to_opt_param_list(opt_dim,params)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE

  integer		:: opt_dim

  real*8		:: params(1:opt_dim) 
 
! local variables

  integer		:: i
  integer		:: param_count
  integer		:: pole_count
  integer		:: zero_count

! START 

! check...
  if (opt_dim.NE.(order*2+1)) then
    write(*,*)'Error in impedance_S_PZ_filter_to_opt_param_list'
    write(*,*)'opt_dim discrepancy',opt_dim,order*2+1
    STOP
  end if
  
  param_count=0

! Gain term
  param_count=param_count+1
  params(param_count)=filter_S_PZ(1)%G
  
! Poles  

  pole_count=0
  
! real poles first

  do i=1,filter_S_PZ(1)%n_real_poles
    param_count=param_count+1
    pole_count=pole_count+1
    params(param_count)=dble(filter_S_PZ(1)%poles(pole_count))
  end do
  
! complex poles second

  do i=1,filter_S_PZ(1)%n_complex_pole_pairs
    param_count=param_count+1
    pole_count=pole_count+1
    params(param_count)=dble(filter_S_PZ(1)%poles(pole_count))
 ! conjugate pole    
    param_count=param_count+1
    pole_count=pole_count+1
    params(param_count)=imag(filter_S_PZ(1)%poles(pole_count))
  end do
  
! Zeros  

  zero_count=0
  
! real zeros first

  do i=1,filter_S_PZ(1)%n_real_zeros
    param_count=param_count+1
    zero_count=zero_count+1
    params(param_count)=dble(filter_S_PZ(1)%zeros(zero_count))
  end do
  
! complex zeros second

  do i=1,filter_S_PZ(1)%n_complex_zero_pairs
    param_count=param_count+1
    zero_count=zero_count+1
    params(param_count)=dble(filter_S_PZ(1)%zeros(zero_count))
 ! conjugate zero    
    param_count=param_count+1
    zero_count=zero_count+1
    params(param_count)=imag(filter_S_PZ(1)%zeros(zero_count))
  end do
    
! check we have set the expected number of parameters
  if (param_count.NE.opt_dim) then
    write(*,*)'Error in impedance_S_PZ_filter_to_opt_param_list'
    write(*,*)'opt_dim discrepancy',opt_dim,param_count
    STOP 
  end if
       
END SUBROUTINE impedance_S_PZ_filter_to_opt_param_list
!
! NAME
!     opt_param_list_to_impedance_S_PZ_filter
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 7/1/2013 CJS
!
!
SUBROUTINE opt_param_list_to_impedance_S_PZ_filter(opt_dim,params)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE

  integer		:: opt_dim

  real*8		:: params(1:opt_dim) 
 
! local variables

  integer		:: i
  integer		:: param_count
  integer		:: pole_count
  integer		:: zero_count

! START 

! check...
  if (opt_dim.NE.(order*2+1)) then
    write(*,*)'Error in opt_param_list_to_impedance_S_PZ_filter'
    write(*,*)'opt_dim discrepancy',opt_dim,order*2+1
    STOP
  end if
  
  param_count=0

! Gain term
  param_count=param_count+1
  filter_S_PZ(1)%G=params(param_count)
  
! Poles  

  pole_count=0
  
! real poles first

  do i=1,filter_S_PZ(1)%n_real_poles
    param_count=param_count+1
    pole_count=pole_count+1
    filter_S_PZ(1)%poles(pole_count)=params(param_count)
  end do
  
! complex poles second

  do i=1,filter_S_PZ(1)%n_complex_pole_pairs
    param_count=param_count+1
    pole_count=pole_count+1
    filter_S_PZ(1)%poles(pole_count)  =DCMPLX(params(param_count), params(param_count+1))
    filter_S_PZ(1)%poles(pole_count+1)=DCMPLX(params(param_count),-params(param_count+1))
    param_count=param_count+1
    pole_count=pole_count+1
  end do
  
! Zeros  

  zero_count=0
  
! real zeros first

  do i=1,filter_S_PZ(1)%n_real_zeros
    param_count=param_count+1
    zero_count=zero_count+1
    filter_S_PZ(1)%zeros(zero_count)=params(param_count)
  end do
  
! complex zeros second

  do i=1,filter_S_PZ(1)%n_complex_zero_pairs
    param_count=param_count+1
    zero_count=zero_count+1
    filter_S_PZ(1)%zeros(zero_count)  =DCMPLX(params(param_count), params(param_count+1))
    filter_S_PZ(1)%zeros(zero_count+1)=DCMPLX(params(param_count),-params(param_count+1))
    param_count=param_count+1
    zero_count=zero_count+1
  end do
  
! check we have set the expected number of parameters
  if (param_count.NE.opt_dim) then
    write(*,*)'Error in opt_param_list_to_impedance_S_PZ_filter'
    write(*,*)'opt_dim discrepancy',opt_dim,param_count
    STOP 
  end if
       
END SUBROUTINE opt_param_list_to_impedance_S_PZ_filter
