       subroutine opt_simplex()
!
! Optimise using the downhill simplex method
!
       
USE opt_module         

IMPLICIT NONE   
       
  real*8 min_error
!  
  real*8,allocatable	   :: point_set(:,:)
  real*8,allocatable	   :: error(:)
  real*8,allocatable	   :: cg(:)
!  
  integer n_calc_points
  integer i,k,id,p
  integer n,m,trialpt,ir,ic,first_pos
!
  integer nhigh,nlow,nshigh,neval,ns,np,ncalc
  real*8 fhigh,flow,ftry,fshigh,fsave
  real*8 cvg,fcvg,tcvg,swap
  real*8 oldfhigh
  integer iswap
  
  real*8 new_point
       
! START

  write(*,*)'Optimisation using the simplex method'
  
! number of calculation points in optimisation calculations    

  fcvg=convergence_criterion
  tcvg=convergence_criterion
  np=n_params
  ns=np+1
  neval=0
 
  n_calc_points=1+n_params
       
! number of calculation points in optimisation calculations       
  
  allocate( point_set(n_params,n_calc_points) )
  allocate( error(n_calc_points) )
  allocate( cg(n_calc_points) )
    
  min_error=1e30

! initial point is at the centre of the point set.
  do p=1,n_params
    point_set(p,1)=parameters(p)
  end do
       
! Set all other points to be the same as point 1 initially   
  do p=1,n_params
    do i=2,n_calc_points
      point_set(p,i)=point_set(p,1)
    end do
  end do
      
! add offset to other points...       
  do p=1,n_params
    CALL get_new_point(point_set(p,1),delta(p),point_set(p,p+1))
  end do

! Evaluate the mean square error at each of the trial points

  do trialpt=1,n_params+1
       
!   ******** Evaluate all points ********
    do p=1,n_params
      parameters(p)=point_set(p,trialpt)
    end do
    call opt_calc_error()
    error(trialpt)=objective_function(evaluation_number)
    if (evaluation_number.ge.max_evaluations) GOTO 1000  
       	 	      
  end do
!
! find highest, second highest and lowest indexes in the simplex
!
  fhigh=1d99
1   continue
    oldfhigh=fhigh
    nhigh=0
    nshigh=0
    nlow=0
    fhigh=0d0
    flow=1d30
    fshigh=0d0
    do 10 i=1,ns
! 
      if(error(i).lt.flow) then
        nlow=i
        flow=error(i)
      end if
      if(error(i).ge.fhigh) then
        nshigh=nhigh
        nhigh=i
        fshigh=fhigh
        fhigh=error(i)
      else if(error(i).ge.fshigh) then
        nshigh=i
        fshigh=error(i)
      end if

10  continue

!
! test for convergence
!
    cvg=2d0*abs((error(nhigh)-error(nlow))/(error(nhigh)+error(nlow)))

    if ((cvg.lt.fcvg).or.(error(nlow).lt.fcvg)) then

! put best solution at simplex origin
      do 20 i=1,np
        swap=point_set(i,1)
        point_set(i,1)=point_set(i,nlow)
        point_set(i,nlow)=swap
20    continue

      swap=error(1)
      error(1)=error(nlow)
      error(nlow)=swap

!      print*,'convergence of simplex procedure,',   &
!            ' number of iterations=',neval
!      print*,'fhigh=',fhigh,' flow=',flow
!      print*,'cvg=',cvg

!   ******** Evaluate at new point ********
      trialpt=1
      do p=1,n_params
        parameters(p)=point_set(p,trialpt)
      end do
      call opt_calc_error()
      error(trialpt)=objective_function(evaluation_number)
      if (evaluation_number.ge.max_evaluations) GOTO 1000  
			
! optimum point is at the centre of the point set.
      min_error=error(1)	
      return

    end if ! converged solution

! calculate centre of gravity of simplex

    do 5 k=1,np
      cg(k)=0d0
        ncalc=0
    	do 6 i=1,ns
    	  if (i.ne.nhigh) then
            ncalc=ncalc+1
    	    cg(k)=cg(k)+point_set(k,i)
    	  end if
6   	continue
      cg(k)=cg(k)/dble(ncalc)
5   continue

! Reflection step; reverse offset of highest point wrt lowest point

    do 80 i=1,np
      delta(i)=point_set(i,nhigh)-cg(i)
      
      CALL get_new_point(cg(i),delta(i),point_set(i,nhigh))
      
80  continue

! Reconfigure simplex

    fhigh=error(nhigh)
    neval=neval+1

!   ******** Evaluate at new point ********
    trialpt=nhigh
    do p=1,n_params
      parameters(p)=point_set(p,trialpt)
    end do
    call opt_calc_error()
    error(trialpt)=objective_function(evaluation_number)
    if (evaluation_number.ge.max_evaluations) GOTO 1000  

    ftry=error(nhigh)

    if (error(nhigh).lt.error(nlow)) then

! if new point is better than the best so far, try a stretch 

      do 90 i=1,np
	
        CALL get_new_point(cg(i),1.9d0*delta(i),point_set(i,nhigh))
	
90    continue
      neval=neval+1

!   ******** Evaluate at new point ********
      trialpt=nhigh
      do p=1,n_params
        parameters(p)=point_set(p,trialpt)
      end do
      call opt_calc_error()
      error(trialpt)=objective_function(evaluation_number)
      if (evaluation_number.ge.max_evaluations) GOTO 1000  

! if this is worse than the original reflection then revert back
      if (error(nhigh).gt.ftry) then
!        write(*,*)'Reflect chosen, stable model'
        error(nhigh)=ftry
        do 95 i=1,np
          CALL get_new_point(cg(i),0.9d0*delta(i),point_set(i,nhigh))
95      continue
      else 
!        write(*,*)'Stretch point chosen'
      end if
!
! else if the reflected point is worse than the second highest,
! try a 1 dimensional contraction
!
    else if (ftry.ge.error(nshigh)) then
      fsave=ftry
      do 100 i=1,np
        CALL get_new_point(cg(i),0.475d0*delta(i),point_set(i,nhigh))
100   continue
      neval=neval+1

!   ******** Evaluate at new point ********
      trialpt=nhigh
      do p=1,n_params
        parameters(p)=point_set(p,trialpt)
      end do
      call opt_calc_error()
      error(trialpt)=objective_function(evaluation_number)
      if (evaluation_number.ge.max_evaluations) GOTO 1000  

! check for improvement, if there is none then do a full contraction
      if (error(nhigh).ge.fsave) then
! restore original simplex
        do 110 i=1,np
          CALL get_new_point(cg(i),delta(i),point_set(i,nhigh))
110     continue

! contraction procedure around minimum value point

!        write(*,*)'ND contraction chosen'
        do 40 i=1,ns
          if (i.ne.nlow) then
            do 50 id=1,np
              delta(id)=point_set(id,i)-point_set(id,nlow)
              CALL get_new_point(point_set(id,nlow),delta(id)/2.1d0,point_set(id,i))
50          continue
            neval=neval+1
	      
!   ******** Evaluate at new point ********
            trialpt=i
            do p=1,n_params
              parameters(p)=point_set(p,trialpt)
            end do
            call opt_calc_error()
            error(trialpt)=objective_function(evaluation_number)
            if (evaluation_number.ge.max_evaluations) GOTO 1000  
          end if
40      continue
      else 
!        write(*,*)'1D contraction chosen'
      end if
    else
!       write(*,*)'Reflect chosen'
    end if

!
    goto 1
!
    return
      
1000 CONTINUE      

  deallocate( point_set )
  deallocate( error )
  deallocate( cg )

  return
      
      end
!      
!      
!
  SUBROUTINE get_new_point(cg,delta,new_point)

IMPLICIT NONE   
  
real*8  cg,delta,new_point
real*8  new_point_t

! START

  new_point_t=cg-delta
  if ( new_point_t.ge.1d0 ) then 
    new_point=0.9d0+0.1d0*cg
  else if (new_point_t.le.0d0) then
    new_point=0.1d0*cg
  else 
    new_point=new_point_t
  end if

  if ( (new_point.lt.0d0).OR.(new_point.gt.1d0) ) then
    write(*,*)'Error: simplex point out of range'
    stop
  end if

  return
  
  END SUBROUTINE get_new_point
