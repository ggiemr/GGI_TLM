! MODULE opt_module

!
! NAME
!     MODULE opt_module
!
! DESCRIPTION
!     data relating to optimisation
!	
!     
! COMMENTS
!     
!
!
!
! HISTORY
!
!     started 19/01/10 CJS
!
!

MODULE opt_module

character		:: optimisation_technique

integer 		:: n_params
real*8,allocatable 	:: parameters(:)
real*8,allocatable 	:: delta(:)

integer 		:: max_evaluations
real*8,allocatable 	:: objective_function(:)
real*8,allocatable 	:: parameter_list(:,:)

integer 		:: evaluation_number
integer 		:: min_evaluation_number
real*8 			:: min_function

integer			:: converged_flag

real*8 			:: convergence_criterion

integer			:: last_call
integer			:: write_values

END MODULE opt_module
