       SUBROUTINE opt_bisect()
!
! Optimise using bisection like routine (Acton).
!
       
USE opt_module            

IMPLICIT NONE   

       real min_error
              
       real*8,allocatable	:: point_set(:,:)
       real*8,allocatable	:: error(:)
       
       integer n_calc_points
       integer i,p,point,iteration
       
       real*8 minimum_error
       real*8 maximum_error,maximum_error_cvg
       integer minimum_error_point
       real*8 last_error,cvg
       integer pmin,pmax
       integer npoles
       real*8 large
       
       real*8 tf,tfp,tfm

! START

  write(*,*)'Optimisation using the bisection method'
  
! number of calculation points in optimisation calculations       
  n_calc_points=1+2*n_params
  
  allocate( point_set(n_params,n_calc_points) )
  allocate( error(n_calc_points) )
 
  large=1e30
                     
! calculate initial set of points to calculate over
       
! initial point is at the centre of the point set.
  do p=1,n_params
    point_set(p,1)=parameters(p)
  end do

  last_error=1e30

! start optimisation loop
  do iteration=1,max_evaluations
       
! Set all other points to be the same as point 1 initially   
    do p=1,n_params
      do i=2,n_calc_points
    	point_set(p,i)=point_set(p,1)
      end do
    end do
      
! add offset to other points...       
    do p=1,n_params
      point_set(p,p*2  )=point_set(p,1)+delta(p)
      point_set(p,p*2+1)=point_set(p,1)-delta(p)
    end do

! calculate RMS error for each of the points       
    do point=1,n_calc_points
    
      do p=1,n_params
        parameters(p)=point_set(p,point)
      end do
      call opt_calc_error()
      error(point)=objective_function(evaluation_number)
      
      if (evaluation_number.ge.max_evaluations) then
! we have done enough evaluations
        goto 1000
      end if
      
    end do
	 
! look at each direction and see whether we should increase or decrease delta
    do p=1,n_params
      tf=error(1)
      tfp=error(2*p)
      tfm=error(2*p+1)     
      if ((tfm.gt.tf).and.(tf.lt.tfp)) then
! we do bracket a minimum so decrease delta in this direction		  
    	delta(p)=delta(p)/2.0
      else
    	delta(p)=delta(p)*1.5
	
        if (point_set(p,p*2  )+delta(p).ge.1d0) then
          delta(p)=0.9d0*(1d0-point_set(p,p*2))
        end if
        if (point_set(p,p*2+1)-delta(p).le.0d0) then
          delta(p)=0.9d0*(point_set(p,p*2+1))
       end if
	
	
      end if
    end do

! work out maximum and minimum error point
    
    minimum_error_point=0
    minimum_error=1e30
    maximum_error=0.0
    maximum_error_cvg=0.0
    do point=1,n_calc_points
      if(error(point).lt.minimum_error) then
    	minimum_error=error(point)
    	minimum_error_point=point
      end if
      if(error(point).gt.maximum_error) then
    	maximum_error=error(point)
    	if (error(point).ne.large)then
    	  maximum_error_cvg=error(point)
    	end if 
      end if
    end do
       
    if(minimum_error_point.ne.1) then
! move point 1 to minimum point	 
      do p=1,n_params
    	point_set(p,1)=point_set(p,minimum_error_point)
      end do
    end if

! return if we are at the minimum       
    cvg=abs(maximum_error_cvg-minimum_error)
    if(cvg.lt.convergence_criterion) then
      converged_flag=1
      goto 1000
    end if
    	 
    if (last_error.ne.minimum_error) then	
      last_error=minimum_error
    end if

! next optimisation iteration
  end do

1000 continue  

  deallocate( point_set )
  deallocate( error )

  return
  
  END
