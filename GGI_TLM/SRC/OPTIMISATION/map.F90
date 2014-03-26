       SUBROUTINE map()
!
! Search optimisation space in a systematic way to create a map
! of the response function
!
       
USE opt_module            

IMPLICIT NONE   

! local variables

  integer :: n_d,n,n_tot
  real*8  :: n_dr
  
  real*8,allocatable :: x(:)
  integer,allocatable :: nx(:)
  real*8             :: dx
  
  integer p

! START

! number of function evaluations in each dimension of the optimisation space
  n_dr=dble(max_evaluations)**(1d0/dble(n_params))
  n_d=int(n_dr)
  n_tot=n_d**n_params
  
  write(*,*)'max_evaluations=',max_evaluations
  write(*,*)'n_dr =',n_dr
  write(*,*)'n_d  =',n_d
  write(*,*)'n_tot=',n_tot
  
  allocate( x(1:n_params) )
  allocate( nx(1:n_params) )
  
  dx=1d0/dble(n_d)
  
  nx(:)=0
  
  write(*,*)'Optimisation using the mapping method'
  
  do n=1,n_tot
    if (n.ne.1) then
      nx(1)=nx(1)+1
    end if
    do p=1,n_params
      if (nx(p).eq.n_d) then
        nx(p)=0
	if (p+1.le.n_params) then
	  nx(p+1)=nx(p+1)+1
	else
	  write(*,*)'Error in map'
	  write(*,*)'n=',n
	  stop
	end if
      end if
    end do
    
    do p=1,n_params
      x(p)=nx(p)*dx+dx/2d0
      parameters(p)=x(p)
    end do
    call opt_calc_error()
      
    if (evaluation_number.ge.max_evaluations) then
! we have done enough evaluations
      goto 1000
    end if
      
  end do
 
1000 continue  

  write(*,*)'max_evaluations=',max_evaluations
  write(*,*)'n_dr =',n_dr
  write(*,*)'n_d  =',n_d
  write(*,*)'n_tot=',n_tot

  deallocate( x )
  deallocate( nx )

  return
  
  END
