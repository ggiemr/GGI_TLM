!
! Module to Evaluate an expression in a character string which may include variables
! Variables are specified by $N where N is an integer
! Example expressions would be '1+3/2' '$1+$2*2+5' '1/($2+1)'
!
!

MODULE evaluate_string_expression

IMPLICIT NONE

! operator types

integer,parameter :: plus=1
integer,parameter :: minus=2
integer,parameter :: times=3
integer,parameter :: divide=4

! variables

integer :: n_vars
real*8,allocatable :: var_list(:)

CONTAINS

RECURSIVE SUBROUTINE eval_expression(expr,value,level_in,verbose)

IMPLICIT NONE

! variables passed to subroutine

character*256,intent(IN) :: expr
real*8,intent(OUT)       :: value
integer,intent(IN)       :: level_in
logical,intent(INOUT)       :: verbose

! local variables

integer :: pos
integer :: length

character*256 :: left_string
character*256 :: right_string
character*256 :: rem_string

integer       :: level
character*256 :: level_string

real*8  :: left_value
real*8  :: right_value

logical :: operator_found
integer           :: operator1
integer           :: operator2

integer :: i

! START

level=level_in+1
level_string=''
do i=1,(level-1)*6
  level_string(i:i)=' '
end do
level_string(i:i)=':'

if(verbose) then
  write(*,*)trim(level_string),'CALLED eval_expression with string:'
  write(*,*)trim(level_string),trim(expr)
  write(*,*)trim(level_string),'level        =',level
end if

length=LEN(trim(expr))

if(verbose) write(*,*)trim(level_string),'string length=',length

if(length.EQ.0) then
  if(verbose) write(*,*)'ERROR: called eval_expression with empty string; returning zero'
  value=0d0
  RETURN
end if

pos=0
operator_found=.FALSE.

! find the next + or - operator and split the string at this point into left and right hand strings
! note that we ignore anything in brackets here

CALL get_next_pm_operator(expr,operator_found,operator1,left_string,right_string,pos,level_string,verbose)

if (operator_found) then
  if (pos.EQ.length) then
    write(*,*)'ERROR: found operator with no right hand value'
    STOP
  end if
  
  if (verbose) then
    write(*,*)trim(level_string),'left string =',trim(left_string)
    write(*,*)trim(level_string),'right string=',trim(right_string)
    write(*,*)trim(level_string),'Evaluating left string...'
  end if
    
! evaluate the left string

  CALL eval_expression(left_string,left_value,level,verbose)
  
  if (verbose) write(*,*)trim(level_string),'left_value:',left_value

! evaluate the plus/ minus terms from left to right  
  do while (operator_found)
  
    rem_string=right_string
    
! extract the next value and operator
  
    if (verbose) then
      write(*,*)trim(level_string),'Extracting the next right hand term from string:'
      write(*,*)trim(level_string),trim(rem_string)
    end if
    CALL get_next_pm_operator(rem_string,operator_found,operator2,left_string,right_string,pos,level_string,verbose)

    if (operator_found) then
    
      if (verbose) then
        write(*,*)trim(level_string),'Found next right hand term as string:'
        write(*,*)trim(level_string),trim(left_string)
      end if

      if (verbose) write(*,*)trim(level_string),'Evaluating next right hand term...'
      
      CALL eval_expression(left_string,right_value,level,verbose)
  
      if (verbose) write(*,*)trim(level_string),'right_value:',right_value
      
    else
    
      if (verbose) then
        write(*,*)trim(level_string),'Last term found as string:'
        write(*,*)trim(level_string),trim(rem_string)
        write(*,*)trim(level_string),'Evaluating last term...'
      end if
      
      CALL eval_expression(rem_string,right_value,level,verbose)
  
      if (verbose) write(*,*)trim(level_string),'right_value:',right_value
    
    end if

! Evaluate using the previously found operator    
    if (operator1.EQ.minus) then
    
      value=left_value-right_value
      
    else 
      
      value=left_value+right_value
      
    end if  ! operator_type
    
    if(operator_found) then
! shift variables for the next operator to be applied
      operator1=operator2
      left_value=value    
    end if
  
  end do ! while operator_found
    
else
  
! no + or - operator has been found so the string should be evaluated from left to right

  if (verbose) write(*,*)trim(level_string),'No + or - operator found. Evaluating string expression:'
  if (verbose) write(*,*)trim(level_string),trim(expr)

! get first value and operator (if it exists)
  CALL get_next_td_operator(expr,operator_found,operator1,left_string,right_string,pos,level_string,verbose)
  CALL get_string_value(left_string,left_value,level,level_string,verbose)
  
  if (verbose) write(*,*)trim(level_string),'left_value:',left_value
  
  if (.NOT.operator_found) then
    value=left_value      
    if (verbose) write(*,*)trim(level_string),'Value=',value
    RETURN
  end if
  
  do while (operator_found)
  
    rem_string=right_string
    
! extract the next value and operator
    CALL get_next_td_operator(rem_string,operator_found,operator2,left_string,right_string,pos,level_string,verbose)
    CALL get_string_value(left_string,right_value,level,level_string,verbose)
  
    if (verbose) write(*,*)trim(level_string),'right_value:',right_value

! Evaluate using the previously found operator    
    if (operator1.EQ.times) then
    
      value=left_value*right_value
      
    else if (operator1.EQ.divide) then
    
      if (right_value.EQ.0d0) then
        write(*,*)'ERROR: division by zero in eval_expression'
        STOP
      end if 
      
      value=left_value/right_value
      
    end if  ! operator_type
    
    if(operator_found) then
! shift variables for the next operator to be applied
      operator1=operator2
      left_value=value    
    end if
    
  end do ! while operator_found
  
  if (verbose) write(*,*)trim(level_string),'Value=',value
  
end if

RETURN

END SUBROUTINE eval_expression
!
! ____________________________________________________
!
!
RECURSIVE SUBROUTINE get_string_value(number,value,level_in,level_string,verbose)

IMPLICIT NONE

! variables passed to subroutine

character*256,intent(IN)      :: number
real*8,intent(OUT)            :: value
integer,intent(IN)            :: level_in
character*256,intent(INOUT)   :: level_string
logical,intent(INOUT)         :: verbose

! local variables

integer :: length
character :: ch
character :: ch2
character*255 :: variable_number
integer :: var

character*256 :: inner_string
integer :: level

! START

level=level_in+1
length=LEN(trim(number))
ch=number(1:1)

if(verbose) then 
  write(*,*)trim(level_string),'CALLED get_string_value with string:'
  write(*,*)trim(level_string),trim(number)
end if

if (ch.EQ.'(') then

! remove the outer brackets from the expresion and evaluate the inner string

! First check that we have a closing bracket at the end of the string
  ch2=number(length:length)
  if (ch2.NE.')') then
    write(*,*)'ERROR: no closing bracket found in expression:',trim(number)
    STOP
  end if
  
  inner_string=''
  inner_string(1:length-2)=number(2:length-1)

! Evaluate the expression in the bracket
  CALL eval_expression(inner_string,value,level,verbose)

else if (ch.EQ.'$') then

  variable_number(1:255)=number(2:256)
  read(variable_number,*,ERR=9010)var
  
  if ( (var.GT.0).AND.(var.LE.n_vars) ) then
    value=var_list(var)
  else
    GOTO 9020
  end if
  
else
! read the numeric value from the string

  read(number,*,ERR=9000)value

end if

RETURN

9000 write(*,*)'ERROR reading value from string:',trim(number)
STOP

9010 write(*,*)'ERROR reading variable number from string:',trim(number)
STOP

9020 write(*,*)'ERROR variable number',var,' has not been set'
STOP

END SUBROUTINE get_string_value
!
! ____________________________________________________
!
!
SUBROUTINE get_next_pm_operator(expr,operator_found,operator,left_string,right_string,pos,level_string,verbose)

! search from the left to find the next plus or minus arithmetic operator where we can split the string +,-

IMPLICIT NONE

! variables passed to subroutine

character*256,intent(IN)  :: expr
logical,intent(OUT)       :: operator_found
integer,intent(OUT)       :: operator
character*256,intent(OUT) :: left_string
character*256,intent(OUT) :: right_string
integer,intent(OUT)       :: pos
character*256,intent(IN) :: level_string
logical,intent(IN) :: verbose

! local variables

integer :: length
integer :: right_string_length
integer :: left_pos,br_count

! START

length=LEN(trim(expr))
pos=0
operator_found=.FALSE.

do while ( (.NOT.operator_found).AND.(pos.LT.length) )
  
  pos=pos+1

  if (expr(pos:pos).EQ.'+') then
  
    operator=plus
    operator_found=.TRUE.
    if(verbose) write(*,*)trim(level_string),'operator + at position ',pos
    
  else if (expr(pos:pos).EQ.'-') then
  
    operator=minus
    operator_found=.TRUE.
    if(verbose) write(*,*)trim(level_string),'operator - at position ',pos
    
  else if (expr(pos:pos).EQ.'(') then
  
    if(verbose) write(*,*)trim(level_string),'opening bracket=( at position ',pos
    
! now read until a closing bracket is found i.e. ignore anything in brackets here
    left_pos=pos
    br_count=1

    do pos=left_pos+1,length
  
      if (expr(pos:pos).EQ.'(') then      
        br_count=br_count+1
      end if
    
      if (expr(pos:pos).EQ.')') then
      
        br_count=br_count-1
      
        if (br_count.Eq.0) then
          if(verbose) write(*,*)'found closing bracket=) at position ',pos        
          exit
        end if
      
      end if
  
    end do
    
  end if    ! character check

end do

if (operator_found) then

  left_string=''
  left_string(1:pos-1)=expr(1:pos-1)
  
  right_string=''
  right_string_length=length-pos
  right_string(1:right_string_length)=expr(pos+1:length)
  
else

  left_string=''
  right_string=''
  
end if

END SUBROUTINE get_next_pm_operator
!
! ____________________________________________________
!
!
SUBROUTINE get_next_td_operator(expr,operator_found,operator,left_string,right_string,pos,level_string,verbose)

! search from the left to find the next times or divide arithmetic operator

IMPLICIT NONE

! variables passed to subroutine

character*256,intent(IN)  :: expr
logical,intent(OUT)       :: operator_found
integer,intent(OUT)       :: operator
character*256,intent(OUT) :: left_string
character*256,intent(OUT) :: right_string
integer,intent(OUT)       :: pos
character*256,intent(IN) :: level_string
logical,intent(IN) :: verbose

! local variables

integer :: length
integer :: right_string_length
integer :: left_pos,br_count

! START

length=LEN(trim(expr))
pos=0
operator_found=.FALSE.

do while ( (.NOT.operator_found).AND.(pos.LT.length) )
  
  pos=pos+1

  if (expr(pos:pos).EQ.'*') then
  
    operator=times
    operator_found=.TRUE.
    if(verbose) write(*,*)trim(level_string),'operator * at position ',pos
    
  else if (expr(pos:pos).EQ.'/') then
  
    operator=divide
    operator_found=.TRUE.
    if(verbose) write(*,*)trim(level_string),'operator / at position ',pos
    
  else if (expr(pos:pos).EQ.'(') then
  
    if(verbose) write(*,*)trim(level_string),'opening bracket=( at position ',pos
    
! now read until a closing bracket is found i.e. ignore anything in brackets here
    left_pos=pos
    br_count=1

    do pos=left_pos+1,length
  
      if (expr(pos:pos).EQ.'(') then      
        br_count=br_count+1
      end if
    
      if (expr(pos:pos).EQ.')') then
      
        br_count=br_count-1
      
        if (br_count.Eq.0) then
          if(verbose) write(*,*)'found closing bracket=) at position ',pos        
          exit
        end if
      
      end if
  
    end do
           
  end if

end do

if (operator_found) then

  left_string=''
  left_string(1:pos-1)=expr(1:pos-1)
  
  right_string=''
  right_string_length=length-pos
  right_string(1:right_string_length)=expr(pos+1:length)
  
else

  left_string=expr
  right_string=''
  
end if

END SUBROUTINE get_next_td_operator
!
! ____________________________________________________
!
!
SUBROUTINE get_rem_string(right_string,operator_found,operator,rem_string,level_string,verbose)

! serch from the left to find the next arithmetic operator, +,- i.e. and split the input string into the right string and the remainder string

IMPLICIT NONE

! variables passed to subroutine

character*256,intent(INOUT):: right_string
logical,intent(OUT)        :: operator_found
integer,intent(OUT)        :: operator
character*256,intent(OUT)  :: rem_string
character*256,intent(IN)   :: level_string
logical,intent(IN)         :: verbose

! local variables

character*256 :: ip_string
integer :: length
integer :: pos
integer :: rem_string_length

! START

ip_string=right_string

length=LEN(trim(ip_string))
pos=0
operator_found=.FALSE.

do while ( (.NOT.operator_found).AND.(pos.LT.length) )
  
  pos=pos+1

  if (ip_string(pos:pos).EQ.'+') then
  
    operator=plus
    operator_found=.TRUE.
    if(verbose) write(*,*)trim(level_string),'operator + at position ',pos
    
  else if (ip_string(pos:pos).EQ.'-') then
  
    operator=minus
    operator_found=.TRUE.
    if(verbose) write(*,*)trim(level_string),'operator - at position ',pos
    
  end if

end do

if (operator_found) then

  right_string=''
  right_string(1:pos-1)=ip_string(1:pos-1)
  
  rem_string=''
  rem_string_length=length-pos
  rem_string(1:rem_string_length)=ip_string(pos+1:length)
  
else

  rem_string=''
  
end if

END SUBROUTINE get_rem_string


END MODULE evaluate_string_expression

