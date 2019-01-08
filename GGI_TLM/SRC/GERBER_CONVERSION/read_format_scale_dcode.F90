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
!     SUBROUTINE read_format
!     SUBROUTINE read_scale
!     SUBROUTINE read_Dcode
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
SUBROUTINE read_format(line,xfni,xfnd,yfni,yfnd)

IMPLICIT NONE

character*256 :: line
integer       :: xfni,xfnd,yfni,yfnd

character*5 :: command
character   :: ch1,ch2
character*2 :: eol

! START

read(line,'(A5,A,I1,I1,A,I1,I1,A2)')command,ch1,xfni,xfnd,ch2,yfni,yfnd,eol

! Check the format is OK
if (command(1:5).NE.'%FSLA') then
  write(*,*)'Expecting:%FSLA',' found:',command(1:5)
  STOP
end if

if (ch1.NE.'X') then
  write(*,*)'Expecting:X',' found:',ch1
  STOP
end if

if (ch2.NE.'Y') then
  write(*,*)'Expecting:Y',' found:',ch2
  STOP
end if

if (eol(1:2).NE.'*%') then
  write(*,*)'Expecting:*%',' found:',eol(1:2)
  STOP
end if

write(*,'(A9,I1,I1,A10,I1,I1)')'X format:',xfni,xfnd,' Y format:',yfni,yfnd

RETURN

END SUBROUTINE read_format
!
! __________________________________________________________________________
!
!
SUBROUTINE read_scale(line,scale)

IMPLICIT NONE

character*256 :: line
real*8        :: scale

character*3 :: command
character*2   :: ch2
character*2 :: eol

! START

read(line,'(A3,A2,A2)')command,ch2,eol

! Check the format is OK
if (command(1:3).NE.'%MO') then
  write(*,*)'Expecting:%MO',' found:',command(1:3)
  STOP
end if

if (ch2(1:2).EQ.'MM') then
  scale=1D-3
  write(*,*)'Scale is mm, s=',scale
else if (ch2(1:2).EQ.'IN') then
  scale=25.4D-3
  write(*,*)'Scale is inches, s=',scale
else
  write(*,*)"Expecting:'MM' or 'IN' "," found:",ch2
  STOP
end if

if (eol(1:2).NE.'*%') then
  write(*,*)'Expecting:*%',' found:',eol(1:2)
  STOP
end if

RETURN

END SUBROUTINE read_scale
!
! __________________________________________________________________________
!
!
SUBROUTINE read_Dcode(line,s,nDcodes,maxDcodes,maxAparams,Dcode,Atype,Aparams,Asize)

IMPLICIT NONE

character*256 :: line
real*8    :: s          ! scale to metres
integer   :: maxDcodes
integer   :: maxAparams
integer   :: nDcodes
integer   :: Dcode(1:maxDcodes)
character :: Atype(1:maxDcodes)
real*8    :: Aparams(1:maxDcodes,1:maxAparams)
real*8    :: Asize(1:maxDcodes)

character*4 :: command
character   :: ch
character*2 :: ch2
character*2 :: eol
integer :: n_params,param
integer :: pos,oldpos

! START

read(line,'(A4)')command

! Check the format is OK
if (command(1:4).NE.'%ADD') then
  write(*,*)'Expecting:%ADD',' found:',command(1:3)
  STOP
end if

nDcodes=nDcodes+1

! read integer Dcode number
pos=5  ! next character position
CALL read_integer_from_line(line,pos,Dcode(nDcodes))

! read template character i.e. the aperture type
Atype(nDcodes)=line(pos:pos)
pos=pos+1

if (Atype(nDcodes).EQ.'C') then
  n_params=2
!  write(*,*)'Dcode=C'
else if (Atype(nDcodes).EQ.'R') then
  n_params=3
!  write(*,*)'Dcode=R'
else if (Atype(nDcodes).EQ.'O') then
  n_params=3
!  write(*,*)'Dcode=O'
else if (Atype(nDcodes).EQ.'P') then
  n_params=4
!  write(*,*)'Dcode=P'
else
  write(*,*)"Expecting one of 'C', 'R', 'O', 'P', found:'",Atype(nDcodes),"'"
  STOP
end if

! read the aperture parameters

Aparams(nDcodes,1:maxAparams)=0d0

if (line(pos:pos).EQ.',') then
pos=pos+1
! read a parameter list

  write(*,*)'Reading parameter list'
  
  do param=1,n_params
  
    oldpos=pos
    
    CALL read_real_from_line(line,pos,Aparams(nDcodes,param))
    
    if (oldpos.EQ.pos) EXIT         ! no number has been read so exit the reading loop
    if (line(pos:pos).NE.'X') EXIT  ! no further parameters to read
    pos=pos+1
    
  end do
end if

if (Atype(nDcodes).EQ.'C') then
  Aparams(nDcodes,:)=Aparams(nDcodes,:)*s
  Asize(nDcodes)= Aparams(nDcodes,1)                        ! circle diameter
else if (Atype(nDcodes).EQ.'R') then
  Aparams(nDcodes,:)=Aparams(nDcodes,:)*s
  Asize(nDcodes)=max(Aparams(nDcodes,1),Aparams(nDcodes,2)) ! max dimension
else if (Atype(nDcodes).EQ.'O') then
  Aparams(nDcodes,:)=Aparams(nDcodes,:)*s
  Asize(nDcodes)=max(Aparams(nDcodes,1),Aparams(nDcodes,2)) ! max dimension
else if (Atype(nDcodes).EQ.'P') then
  Aparams(nDcodes,1)=Aparams(nDcodes,1)*s
  Aparams(nDcodes,4)=Aparams(nDcodes,4)*s
  Asize(nDcodes)= Aparams(nDcodes,1)                        ! circum circle diameter
end if

eol=line(pos:pos+1)
if (eol.NE.'*%') then
  write(*,*)'Expecting:*%',' found:',eol
  STOP
end if

RETURN

END SUBROUTINE read_Dcode
