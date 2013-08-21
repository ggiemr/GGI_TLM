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
! SUBROUTINE optimise_filter
! SUBROUTINE bisection
! SUBROUTINE calc_error
!
! NAME
!     optimise_filter
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
SUBROUTINE optimise_filter

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE filter_functions
USE constants

IMPLICIT NONE
 
! local variables
  
  integer		:: opt_dim

  real*8,allocatable	:: params(:) 
  
  integer		:: function_loop
  
! START

! setup the optimisation

! Transform the filter(s) into Pole-Zero format  

  do function_loop=1,n_functions
    filter_S_PZ(function_loop)=Convert_filter_S_to_S_PZ( filter_S(function_loop) )
  end do ! next function

! number of optimisation parameters

  if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
  
    opt_dim=order*2+2  ! 2*order+1 filter coefficients plus a conductivity

    ALLOCATE( params(1:opt_dim) )
    CALL dielectric_magnetic_S_PZ_filter_to_opt_param_list(opt_dim,params)
       
  else if (fit_type.eq.thin_layer) then  
    
    opt_dim=3*(order*2+1)  ! 3*(2*order+1) filter coefficients (note 3* as Z12=Z21)
    ALLOCATE( params(1:opt_dim) )
    CALL thin_layer_S_PZ_filter_to_opt_param_list(opt_dim,params)
    
  else if (fit_type.eq.impedance) then  
    
    opt_dim=order*2+1  ! 2*order+1 filter coefficients 

    ALLOCATE( params(1:opt_dim) )
    CALL impedance_S_PZ_filter_to_opt_param_list(opt_dim,params)
    
  end if

! CALL the optimisation routine...  
  CALL bisection(opt_dim,params)

  if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
  
    CALL opt_param_list_to_dielectric_magnetic_S_PZ_filter(opt_dim,params)
       
  else if (fit_type.eq.thin_layer) then  
    
    CALL opt_param_list_to_thin_layer_S_PZ_filter(opt_dim,params)
    
  else if (fit_type.eq.impedance) then  
    
    CALL opt_param_list_to_impedance_S_PZ_filter(opt_dim,params)
    
  end if
  
  DEALLOCATE( params )

! Transform the filter back into Rational function format  

  do function_loop=1,n_functions
    filter_S(function_loop)=Convert_filter_S_PZ_to_S( filter_S_PZ(function_loop) )
  end do ! next function
  
! Calculate the mean square error of the Pole-Residue format filter as a check on the transformation
 
  CALL calculate_error_Sfilter()
  CALL write_line_real('Final Stabilised Sfilter.    Mean square error=',Mean_square_error,0,ff_output_to_screen)
       
END SUBROUTINE optimise_filter
!
! NAME
!     bisection
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
SUBROUTINE bisection(opt_dim,params)

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
  
  integer		:: n_calc_points

  real*8,allocatable	:: delta(:) 
  real*8,allocatable	:: point_set(:,:)
  real*8,allocatable	:: error(:)
  logical,allocatable	:: stable_system(:)
  
  real*8	:: min_error
  integer	:: i,p,point,iteration
  
  real*8	:: minimum_error
  real*8	:: maximum_error,maximum_error_cvg
  integer	:: minimum_error_point
  real*8	:: last_error,cvg
  integer	:: pmin,pmax
  
  real*8	:: tf,tfp,tfm

! START

! number of calculation points in optimisation calculations       
  n_calc_points=1+2*opt_dim

  ALLOCATE( delta(opt_dim) )
  ALLOCATE( point_set(opt_dim,n_calc_points) )
  ALLOCATE( error(n_calc_points) )
  ALLOCATE( stable_system(n_calc_points) )
  
  delta(opt_dim)=0d0
  point_set(opt_dim,n_calc_points) =0d0
  error(n_calc_points)=0d0
  stable_system(n_calc_points)=.FALSE.
  
! get initial offsets for optimisation       
  do p=1,opt_dim
    if(params(p).ne.0d0) then
      delta(p)=params(p)/200.0D0
    else
! Zero parameter found
      delta(p)=0D0
    end if
  end do
       
! calculate initial set of points to calculate over
       
! initial point is at the centre of the point set.
  do p=1,opt_dim
    point_set(p,1)=params(p)
  end do

  last_error=1e30

! start optimisation loop
  do iteration=1,max_opt_iterations
       
! Set all other points to be the same as point 1 initially   
    do p=1,opt_dim
      do i=2,n_calc_points
    	point_set(p,i)=point_set(p,1)
      end do
    end do
      
! add offset to other points...       
    do p=1,opt_dim
      point_set(p,p*2  )=point_set(p,1)+delta(p)
      point_set(p,p*2+1)=point_set(p,1)-delta(p)
    end do

! calculate RMS error for each of the points       
    do p=1,n_calc_points
      CALL calc_error(point_set,error,stable_system,opt_dim,n_calc_points,p)
    end do
    
    if (iteration.eq.1) then
      write(*,'(A31,F16.6)')'First Iteration:         error=',error(1)  !,' cvg=',cvg
    else
      write(*,'(A)',advance='no')char(13)
      write(*,'(A10,I14,A7,F16.6)',advance='no')'Iteration ',iteration,' error=',error(1)  !,' cvg=',cvg
    end if
	 
! at lest the origin point should give a stable model...
    if (.NOT.stable_system(1)) then

      write(*,*)'Unstable model at the origin of simplex'
      write(*,*)'Iteration=',iteration,error(1),stable_system(1)
      
      write(*,*)'Parameters:'
      do i=1,opt_dim
    	write(*,*)i,params(i)
      end do   
       
      STOP
      
    end if

! set all unstable models in simplex to have very high error
    do p=1,n_calc_points
      if (.NOT.stable_system(p)) then
    	error(p)=large
      end if  
    end do

! look at each direction and see whether we should increase or decrease delta
    do p=1,opt_dim
      tf=error(1)      ! centre point
      tfp=error(2*p)   ! centre+delta(p)
      tfm=error(2*p+1) ! centre-delta(p)
      if ((tfm.gt.tf).and.(tf.lt.tfp)) then
! we do bracket a minimum so decrease delta in this direction		  
    	delta(p)=delta(p)/2.0D0
      else
! enlarge delta in this direction by a small proportion	  
    	delta(p)=delta(p)*1.2D0
      end if
    end do

! work out maximum and minimum error point
    
    minimum_error_point=0
    minimum_error=1e32
    maximum_error=0.0
    maximum_error_cvg=0.0
    do point=1,n_calc_points
      if(error(point).lt.minimum_error) then
    	minimum_error=error(point)
    	minimum_error_point=point
      end if
      if(error(point).gt.maximum_error) then
    	maximum_error=error(point)
      end if
      if (error(point).ne.large)then
    	if(error(point).gt.maximum_error_cvg) then
    	  maximum_error_cvg=error(point)
    	end if
      end if 
    end do
    
    if(minimum_error_point.ne.1) then
! move point 1 to minimum point	 
      do p=1,opt_dim
    	point_set(p,1)=point_set(p,minimum_error_point)
      end do
    end if

! return if error is close to zero       
    if(abs(minimum_error).lt.opt_accuracy) GOTO 100
    
! return if we are at the minimum       
    cvg=abs(maximum_error_cvg-minimum_error)
    if(cvg.lt.opt_accuracy) GOTO 100
    		 
    if (last_error.ne.minimum_error) then	
      last_error=minimum_error
    end if

! next optimisation iteration
  end do

100 CONTINUE ! Jump here when converged or we have completed all iterations

  write(*,*)
  write(*,'(A31,F16.6,A4,F16.6)')'Final Iteration,         error=',error(1),' cvg=',cvg

! recalculate frequency response and error - this has the effect of
! putting the optimum solution for the zeros into the c array plus the
! conductivity term if it is a material optimisation       	
  p=1 
  CALL calc_error(point_set,error,stable_system,opt_dim,n_calc_points,p)
       
! optimum point is at the centre of the point set.
  do p=1,opt_dim
    params(p)=point_set(p,1)
  end do
  min_error=error(1)	   

  DEALLOCATE( delta )
  DEALLOCATE( point_set )
  DEALLOCATE( error )
  DEALLOCATE( stable_system )
       
END SUBROUTINE bisection
!
! NAME
!     bisection
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
SUBROUTINE calc_error(point_set,error,stable_system,opt_dim,n_calc_points,p)

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff
USE constants

IMPLICIT NONE

  integer		:: opt_dim,n_calc_points,p

  real*8		:: point_set(1:opt_dim,1:n_calc_points)
  real*8		:: error(1:n_calc_points)
  logical		:: stable_system(1:n_calc_points)
  
! local variables

  real*8,allocatable	:: params(:) 
  

! START

  ALLOCATE( params(1:opt_dim) )

  params(1:opt_dim)=point_set(1:opt_dim,p)
  
  if ( (fit_type.eq.dielectric_material).OR.(fit_type.eq.magnetic_material) ) then 
  
    CALL opt_param_list_to_dielectric_magnetic_S_PZ_filter(opt_dim,params)
       
  else if (fit_type.eq.thin_layer) then  
    
    CALL opt_param_list_to_thin_layer_S_PZ_filter(opt_dim,params)
    
  else if (fit_type.eq.impedance) then  
    
    CALL opt_param_list_to_impedance_S_PZ_filter(opt_dim,params)
    
  end if

  CALL calculate_error_Sfilter_PZ()
  error(p)=Mean_square_error
  
  CALL test_stability_PZ()
  stable_system(p)=stable_filter

  DEALLOCATE( params )
       
END SUBROUTINE calc_error
