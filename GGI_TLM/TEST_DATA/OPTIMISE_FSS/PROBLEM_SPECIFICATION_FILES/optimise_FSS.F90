PROGRAM optimise_FSS

USE opt_module      

IMPLICIT NONE

! START

  open(unit=10,file='progress')
  write(10,*)'STARTED: optimise_FSS'
  close(unit=10)

  write(*,*)'Optimisation of FSS' 

  CALL optimise()
 
  open(unit=10,file='progress')
  write(10,*)'FINISHED: optimise_FSS'
  close(unit=10)
 

END PROGRAM optimise_FSS
!
!
!
SUBROUTINE calculate_function(value)

! evaluate the function to minimise using the optimisation procedure

! 1. Turn the normalised parameters into FSS dimensions, 
! 2. Run and post process the TLM result
! 3. Calculate the objective function value from the post processed data

USE opt_module      
USE Constants     

IMPLICIT NONE

real*8 	:: value

real*8 	:: cell_size,dim1,dim2,dim3,dim4
real*8 	:: t1,t2,t3,t4

real*8	:: r(9)
integer :: i,ii
integer :: nvalues

! START
    
! set the FSS dimensions from the parameter values

  t1=4e-2*parameters(1)
  t2=1e-2*parameters(2)
  t3=2e-2*parameters(3)
  t4=2e-2*parameters(4)

  dim4=t1
  dim3=t1+t2
  dim2=t1+t2+t3
  dim1=t1+t2+t3+t4
  cell_size=t1+t2+t3+t4+2e-2*parameters(5)

! evaluate the value to minimise

! First write the dimensions to a file
  open(unit=50,file='FSS_dimensions')
  write(50,8000)cell_size,dim1,dim2,dim3,dim4
8000 format(5E16.6)
  close(unit=50)
  
! Run the script which runs and post processes the TLM solution
  call system("run_GGI_TLM_FSS_solution run_seq OPTIMISE_FSS")

! Calculate the objective value to optimise from the output result(s)

  value=0d0

  open(unit=50,file='S21.fout')
! read header line
  read(50,*)
! read the data file
  nvalues=0
1000  CONTINUE
      read(50,*,end=1010)(r(ii),ii=1,9)
      if ( (r(1).GE.880E6).AND.(r(1).LE.920E6) ) then
        value=value+r(9)  ! sum the transmission coefficient in dB
	nvalues=nvalues+1
      else if ( (r(1).GE.1780E6).AND.(r(1).LE.1820E6) ) then
        value=value+r(9)  ! sum the transmission coefficient in dB    
	nvalues=nvalues+1
     end if  
     GOTO 1000
    
1010 CONTINUE

     if (nvalues.EQ.0) then
       value=large
     else
       value=value/nvalues
     end if

  close(unit=50)
           
 
END SUBROUTINE calculate_function
