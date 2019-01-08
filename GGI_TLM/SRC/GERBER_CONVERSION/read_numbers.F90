!
!    GGI_TLM Time domain electromagnetic field solver based on the TLM method
!    Copyright (C) 2013  Chris Smartt
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.   
!
!
! NAME
!     SUBROUTINE read_integer_from_line
!     SUBROUTINE read_formatted_coordinate_from_line
!     SUBROUTINE read_real_from_line
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 10/12/19 CJS
!     
!
SUBROUTINE read_integer_from_line(line,pos,inumber)

IMPLICIT NONE

character*256 :: line
integer :: pos,inumber

character*256 :: number
character :: ch
integer :: i,endpos

! START

inumber=0   ! default value

!write(*,*)'Reading integer from:',line(pos:len_trim(line))
!write(*,*)pos,len_trim(line)

do i=pos,len_trim(line)
  
  ch=line(i:i)
!  write(*,*)'ch=',ch
  
  if ( (ichar(ch).LT.ichar('0')) .OR. (ichar(ch).GT.ichar('9')) ) then
! the next character is not a number so exit
    endpos=i-1
    EXIT
  end if    
  endpos=i
  
end do

if (endpos-pos.LT.0) RETURN   ! no number to read

number=line(pos:endpos)
read(number,*)inumber

!write(*,*)'Read integer',trim(number),' ',inumber

! update the line position
pos=endpos+1

RETURN

END SUBROUTINE read_integer_from_line
!
! __________________________________________________________________________
!
!
SUBROUTINE read_formatted_coordinate_from_line(line,pos,x,xfni,xfnd,s)

IMPLICIT NONE

character*256 :: line
integer :: pos
real*8  :: x
integer :: xfni,xfnd
real*8  :: s

integer :: sign

character*256 :: number
character :: ch
integer :: npos,nchar
integer :: dcount
integer :: i

! START

! read optional sign
sign=1
if (line(pos:pos).EQ.'-') then
  sign=-1
  pos=pos+1
else if (line(pos:pos).EQ.'+') then
  sign=1
  pos=pos+1
end if

! read the line character by character building a decimal number string

number=''
nchar=xfni+xfnd

! there may be less than the specified number of decimal digits here...
! if this is the case we pad with leading zeros
npos=0
do i=pos,pos+nchar-1
  
  ch=line(i:i)
    
  if ( (ichar(ch).LT.ichar('0')) .OR. (ichar(ch).GT.ichar('9')) ) then
! the next character is not a number so exit with dcount = number of digits read
    EXIT
  end if    
  
  npos=npos+1
  number(npos:npos)=ch
  
end do

pos=pos+npos

! The length of the integer is npos so pad it out with leading zeros to length nchar

! calculate the number of characters to add
dcount=nchar-npos 

do i=1,npos
! shift right by dcount digits
  number(nchar-i+1:nchar-i+1)=number(npos-i+1:npos-i+1)
end do

! add leading zeros
do i=1,dcount
  number(i:i)='0'
end do

! shift the decimal part to the right one character to make space for the decimal point
nchar=nchar+1
do i=1,xfnd
  number(nchar-i+1:nchar-i+1)=number(nchar-i:nchar-i)
end do

! add the decimal point
number(xfni+1:xfni+1)='.'

read(number,*)x

!write(*,*)'number=',trim(number),' x=',x

x=dble(sign)*x*s

RETURN

END SUBROUTINE read_formatted_coordinate_from_line
!
! __________________________________________________________________________
!
!
SUBROUTINE read_real_from_line(line,pos,rnumber)

IMPLICIT NONE

character*256 :: line
integer :: pos
real*8  :: rnumber

character*256 :: number
character :: ch
integer :: i,endpos

! START

rnumber=0d0  ! default value

!write(*,*)'Reading real from:',line(pos:len_trim(line))
!write(*,*)pos,len_trim(line)

do i=pos,len_trim(line)
  
  ch=line(i:i)
!  write(*,*)'ch=',ch
  
  if (  ( (ichar(ch).LT.ichar('0')) .OR. (ichar(ch).GT.ichar('9')) ).AND.   &
          (ichar(ch).NE.ichar('.')) ) then
! the next character is not a number so exit
    endpos=i-1
    EXIT
  end if    
  endpos=i
  
end do

if (endpos-pos.LT.0) RETURN   ! no number to read

number=line(pos:endpos)  
!write(*,*)'number=',number
read(number,*)rnumber

!write(*,*)'Read real: ',trim(number),' : ',rnumber

! update the line position
pos=endpos+1

RETURN

END SUBROUTINE read_real_from_line
