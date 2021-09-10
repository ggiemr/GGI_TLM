!
! NAME
!     SUBROUTINE optimise
!
! DESCRIPTION
!     general optimisation routines
!	
!     
! COMMENTS
!     
! optimisation procedures to link with other code to calculate
! the function to minimise...
!
! We define the number of normalised (to lie between 0 and 1) parameters
! together with a starting point in the optimisation space and a delta value.
!
! A termination criterion and a maximum number of function calls is set.
!
! The optimsation proceeds by calling the function which calculates the 
! objective function with parameters set by the chosen optimisation technique
! as may times as necessary until the termination criterion is met. 
!       
!
! HISTORY
!
!     started 19/01/10 CJS
!
!

SUBROUTINE  optimise
       
USE opt_module      

IMPLICIT NONE   

! local variables

integer parameter 
integer number_of_evaluations_performed     

! START
       
! read optimisation data

  write(*,*)'Multi-dimensional optimisation'
  write(*,*)''
  write(*,*)'Minimises an objective function with respect to'
  write(*,*)'a number of normalised real parameters'
  write(*,*)'The normalised parameters lie in the range 0<p<1'
  write(*,*)' '
  write(*,*)'The optimisation techniques available are:'
  write(*,*)'bisection'
  write(*,*)'simplex'
  write(*,*)' '
 
  write(*,*)'Enter the number of parameters (maximum of 20):'
  read(*,*)n_params
  write(*,*)'Number of parameters=',n_params
  if (n_params.gt.20) then
    write(*,*)'Maximum number of parameters=20'
    stop
  end if
  
  write(*,*)' '
  write(*,*)'Enter the optimisation technique to apply:'
  write(*,*)'The optimisation techniques available are:'
  write(*,*)'bisection'
  write(*,*)'simplex'
  write(*,*)'map'
  read(*,'(A1)')optimisation_technique
  write(*,*)'Optimisation technique=',optimisation_technique
  
  write(*,*)' '
  write(*,*)'Enter the convergence criterion:'
  read(*,*)convergence_criterion
  write(*,*)'Convergence_criterion=',convergence_criterion
  
  write(*,*)' '
  write(*,*)'Enter the maximum number of function evaluations to use'
  read(*,*)max_evaluations

! if the maximum number of evaluations is set to 1 then assume that we 
! only want a single evaluation for testing, otherwise
! increase the number of evaluations so that we can get at least one 
! round of optimisation after the initial evaluation of the simplex points. 

  if (optimisation_technique.eq.'s') then
  
    if ( (max_evaluations.lt.n_params+1).AND.(max_evaluations.NE.1) ) then
      write(*,*)'Increasing maximum number of evaluations'
      max_evaluations=n_params+1
    end if

  else if  (optimisation_technique.eq.'b') then
  
    if ( (max_evaluations.lt.n_params*2+1).AND.(max_evaluations.NE.1) ) then
      write(*,*)'Increasing maximum number of evaluations'
      max_evaluations=n_params*2+1
    end if

  else if  (optimisation_technique.eq.'m') then

    if ( (max_evaluations.lt.1).AND.(max_evaluations.NE.1) ) then
      write(*,*)'Increasing maximum number of evaluations'
      max_evaluations=1
    end if

  else
  
    write(*,*)'Unknown optmisation technique'
    stop
    
  end if

  write(*,*)'Maximum number of function evaluations=',max_evaluations
  
! allocate memory

  allocate( parameters(1:n_params) )
  allocate( delta(1:n_params) )
  allocate( objective_function(1:max_evaluations) )
  allocate( parameter_list(1:max_evaluations,1:n_params) )
  
  write(*,*)'Enter the intial values for the normalised parameters and variation'
  
  do parameter=1,n_params
  
10  continue

    write(*,8000)'Enter value for parameter ',parameter
8000 format(A26,I4)
    read(*,*)parameters(parameter),delta(parameter)
    
    if ( (parameters(parameter).le.0d0).OR.(parameters(parameter).ge.1d0) ) then
      write(*,*)'Parameter should be in the range 0<p<1'
      goto 10
    end if
    
    if (parameters(parameter)+delta(parameter).ge.1d0) then
      delta(parameter)=0.5d0*(1d0-parameters(parameter))
    end if
    if (parameters(parameter)-delta(parameter).le.0d0) then
      delta(parameter)=0.5d0*(parameters(parameter))
    end if
    
  end do ! next parameter
  
  write(*,*)' '
  write(*,*)'       Parameter             Delta   '
  
  do parameter=1,n_params

    write(*,*)parameters(parameter),delta(parameter)
    
  end do ! next parameter
    
  write(*,*)' '
  
! open monitoring file  
  open(unit=55,file='optimisation_progress')
  
  evaluation_number=0
  min_evaluation_number=0
  min_function=1d30
  converged_flag=0
  write_values=0
  last_call=0
  
  if (optimisation_technique.eq.'s') then
    call opt_simplex() 
  else if  (optimisation_technique.eq.'b') then
    call opt_bisect() 
  else if  (optimisation_technique.eq.'m') then
    call map() 
  else
    write(*,*)'Unknown optmisation technique'
    stop
  end if
  
  number_of_evaluations_performed = evaluation_number
  
  write(*,*)
  write(*,*)
  write(*,*)'Number of function evaluations performed=',number_of_evaluations_performed
  write(*,*)'Optimum point at evaluation number=',  min_evaluation_number
  write(*,*)''
  write(*,*)'Optimum point in normalised optimisation space:'
  
  do parameter=1,n_params

    write(*,*)parameter,parameter_list(min_evaluation_number,parameter)
    
  end do ! next parameter
  
  write(*,*)' '
  write(*,*)'Optimum parameter values'
  
  write_values=1
  
  do parameter=1,n_params
    parameters(parameter)=parameter_list(min_evaluation_number,parameter)
  end do
  evaluation_number=min_evaluation_number-1
  last_call=1

!  if (min_evaluation_number.NE.number_of_evaluations_performed) then
! recalculate the optimum solution so that the optimum solution is the active one.
!    call opt_calc_error()
!  end if

! Always recalculate the function at the optimum point

  call opt_calc_error()
  
! close monitoring file  
  close(unit=55)
  
  write(*,*)''
  write(*,*)'Initial objective function value:',objective_function(1)
  write(*,*)'Optimum objective function value:',objective_function(min_evaluation_number)
  write(*,*)''
  write(*,*)'Number of function evaluations:',number_of_evaluations_performed 
    
! deallocate memory

  deallocate( parameters )
  deallocate( delta )
  deallocate( objective_function )
  deallocate( parameter_list )
             
END
       
