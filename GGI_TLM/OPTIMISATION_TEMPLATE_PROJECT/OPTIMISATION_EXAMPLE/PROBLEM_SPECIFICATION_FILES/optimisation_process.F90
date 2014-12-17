PROGRAM optimisation_process

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
 

END PROGRAM optimisation_process
!
! ______________________________________________________________________
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

character*16 number_string

! START
    
! SET THE PROBLEM PARAMETERS FROM THE NORMALISED PARAMETER VALUES

  t1=4e-2*parameters(1)
  t2=1e-2*parameters(2)
  t3=2e-2*parameters(3)
  t4=2e-2*parameters(4)

  dim4=t1
  dim3=t1+t2
  dim2=t1+t2+t3
  dim1=t1+t2+t3+t4
  cell_size=t1+t2+t3+t4+2e-2*parameters(5)

! Write the problem parameters to the parameter_definition.dat file
! note that it is best not to have any space between the '=' and the number
! due to the manner in which we do the parameter replacement in the input file
! hence the use of the subroutine write_parameter_to_file

  open(unit=50,file='PROBLEM_SPECIFICATION_FILES/parameter_definition.dat')
  
  CALL write_parameter_to_file('#FSS_UNIT_CELL_SIZE=',cell_size,50)
  CALL write_parameter_to_file('#DIM1=',dim1,50)
  CALL write_parameter_to_file('#DIM2=',dim2,50)
  CALL write_parameter_to_file('#DIM3=',dim3,50)
  CALL write_parameter_to_file('#DIM4=',dim4,50)

  close(unit=50)
  
! RUN THE SCRIPT WHICH RUNS AND POST PROCESSES THE GGI_TLM SOLUTION

  call system("run_GGI_TLM_trial_solution")
  
! CALCULATE THE OBJECTIVE VALUE TO OPTIMISE FROM THE OUTPUT RESULT(S)

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
!
! ______________________________________________________________________
!
!
SUBROUTINE write_parameter_to_file(parameter_name_string,value,unit)

! wire the specified parameter name and value to file

IMPLICIT NONE

character*(*) parameter_name_string
real*8 	:: value
integer :: unit

! local variables

character*16 number_string

! START
  
  write (number_string,'(E16.6)')value
  write(unit,8000)trim(parameter_name_string),adjustl(number_string)
8000 format(A,A16)

END SUBROUTINE write_parameter_to_file
