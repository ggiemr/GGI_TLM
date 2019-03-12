SUBROUTINE replace_in_string(line,find_string,replace_string)

IMPLICIT NONE

character(LEN=256),INTENT(INOUT) :: line
character(LEN=*),INTENT(IN)      :: find_string
character(LEN=256),INTENT(IN)    :: replace_string

! local variables

integer :: found_pos
integer :: len_line
integer :: len_found
character(LEN=256) :: left_string
character(LEN=256) :: right_string
character(LEN=256) :: opline

! START

!write(*,*)''
!write(*,*)'CALLED:replace_in_string'
!write(*,*)'LINE   :',trim(line)
!write(*,*)'FIND   :',trim(find_string)
!write(*,*)'REPLACE:',trim(replace_string)

len_line=len(trim(line))
len_found=len(trim(find_string))
found_pos=index(line,trim(find_string))

do while (found_pos.NE.0) 

! replace the found string with the replace string
  left_string=''
  left_string=line(1:found_pos-1)
  right_string=''
  right_string=line(found_pos+len_found:len_line)
  
  opline=''
  opline=trim(left_string)//trim(replace_string)//trim(right_string)
  
  line=''
  line=opline
  
  found_pos=index(line,trim(find_string))

end do

END SUBROUTINE replace_in_string
