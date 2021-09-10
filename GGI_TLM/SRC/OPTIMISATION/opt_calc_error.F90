!
!
! NAME
!      SUBROUTINE opt_calc_error()
!     
!
! DESCRIPTION
!     
!	
!     
! COMMENTS
!
!       
!
! HISTORY
!
!     started 19/01/10 CJS
!
!

SUBROUTINE opt_calc_error()
       
USE opt_module      

IMPLICIT NONE
       
! local variables

  integer param       
  real*8 function
  integer n_params_local,p1,p2,parameter
       
! START

  evaluation_number=evaluation_number+1
  
  call calculate_function(function)
        
  do param=1,n_params
    parameter_list(evaluation_number,param)=parameters(param)
  end do ! next parameter
  objective_function(evaluation_number)=function
  
  write(55,8010)(parameter_list(evaluation_number,parameter),parameter=1,n_params),objective_function(evaluation_number)
8010 format(100E14.4)  
  
  if(function.lt.min_function) then
    min_function=function
    min_evaluation_number=evaluation_number
  end if

  if (write_values.eq.1) return
   
  write(6,'(A)',advance='no')char(13)
 
  n_params_local=n_params
  p1=1
  p2=n_params_local
  if (n_params_local.ge.10) then
    p2=10
    write(6,80010)parameters(p1:p2)
    p1=p2+1
    p2=n_params
    n_params_local=n_params_local-10
  end if
  
  if (n_params_local.eq.1) then
    write(6,8001,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.2) then
    write(6,8002,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.3) then
    write(6,8003,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.4) then
    write(6,8004,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.5) then
    write(6,8005,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.6) then
    write(6,8006,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.7) then
    write(6,8007,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.8) then
    write(6,8008,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.9) then
    write(6,8009,advance='no')parameters(p1:p2),function
  else if (n_params_local.eq.10) then
    write(6,80010,advance='no')parameters(p1:p2),function 
  end if
  
  if (evaluation_number.EQ.1) then
    write(6,*)
  end if
  
  if (last_call.EQ.1) then
    write(6,*)
  end if
  
8001  format(  F8.5,X,E10.4)
8002  format( 2F8.5,X,E10.4)
8003  format( 3F8.5,X,E10.4)
8004  format( 4F8.5,X,E10.4)
8005  format( 5F8.5,X,E10.4)
8006  format( 6F8.5,X,E10.4)
8007  format( 7F8.5,X,E10.4)
8008  format( 8F8.5,X,E10.4)
8009  format( 9F7.4,X,E10.4)
80010 format(10F7.4,X,E9.3)
  
  return
     
END SUBROUTINE opt_calc_error    
