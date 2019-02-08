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
!     SUBROUTINE read_Aperture_Macro
!     SUBROUTINE read_modifiers
!     SUBROUTINE read_variable_definition
!
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 8/2/19 CJS
!     
!

SUBROUTINE read_Aperture_Macro(line)

USE gerber

IMPLICIT NONE

character*256 :: line

! local variables
integer :: length
integer :: pos
character*256 :: name

integer :: i
character :: ch,ch2
integer :: ivar_number

! START

write(*,*)'Reading aperture macro from line:'
write(*,*)line
n_AM=n_AM+1

if (n_AM.GT.maxAM) then
  write(*,*)'ERROR maximum number of aperture macros has been exceeded'
  STOP
end if

length=LEN(trim(line))

! Read aperture macro name

pos=4   ! first character of name
i=0
AMname(n_AM)=''
do pos=4,length
  ch=line(pos:pos)
  if (ch.EQ.'*') exit 
  i=i+1
  AMname(n_AM)(i:i)=line(pos:pos)
end do

if(pos.EQ.length+1) then
  write(*,*)'ERROR reading name of aperture macro:'
  write(*,*)trim(line)
  STOP
end if

write(*,*)'aperture macro name is:',trim(AMname(n_AM))

! read primitives and variable definitions

do while (pos.LT.length)

  pos=pos+1
  ch=line(pos:pos)
      
  if (ch.EQ.'%') then
  
    write(*,*)'End of aperture macro'
    
  else 
  
! read the next primitive data...

    AM_n_primitives(n_AM)=AM_n_primitives(n_AM)+1
    total_AM_primitives=total_AM_primitives+1
    AM_to_primitive_list(n_AM,AM_n_primitives(n_AM))=total_AM_primitives
    AM_primitive_number(total_AM_primitives)=0
    
    if (ch.EQ.'$') then
  
      write(*,*)'Read aperture macro variable definitions'
        
! read the variable number first and set the primitive number to minus the variable number to be defined
      pos=pos+1
      write(*,*)'Reading variable number from position',pos,' line:'
      write(*,*)trim(line)
      CALL read_integer_from_line(line,pos,ivar_number)
      write(*,*)'Reading definition of variable number',ivar_number
      
      AM_primitive_number(total_AM_primitives)=-  ivar_number

! read the variable definition string into the first modifier for this primitive      
      CALL read_variable_definition(line,length,pos)
      write(*,*)'Read variable definition:',trim( AM_modifiers(total_AM_primitives,1) )
  
    else if (ch.EQ.'0') then
  
      write(*,*)'Read comment primitive'
    
! read to the end of the block signified by '*'
      do while (ch(pos:pos).NE.'*')
        pos=pos+1
        if (pos.GT.length) then
          write(*,*)'Error reading comment from aperture macro:',trim(line)
          STOP
        end if
        ch=line(pos:pos)
      end do
    
    else if (ch.EQ.'1') then
  
      write(*,*)'Read circle primitive'
    
      AM_primitive_number(total_AM_primitives)=1    
    
      CALL read_modifiers(line,length,pos)
    
    else if (ch.EQ.'2') then
  
      write(*,*)'Read line primitive'
      pos=pos+1
      ch2=line(pos:pos)
    
      if (ch2.EQ.'0') then
  
        write(*,*)'Read vector line '
    
        AM_primitive_number(total_AM_primitives)=20
    
        CALL read_modifiers(line,length,pos)    
    
      else if (ch2.EQ.'1') then
  
        write(*,*)'Read centre line '
    
        AM_primitive_number(total_AM_primitives)=21
    
        CALL read_modifiers(line,length,pos)    
      
      end if
    
    else if (ch.EQ.'4') then
  
      write(*,*)'Read outline primitive'  
    
      AM_primitive_number(total_AM_primitives)=4
    
      CALL read_modifiers(line,length,pos)
    
    else if (ch.EQ.'5') then
  
      write(*,*)'Read polygon primitive'
    
      AM_primitive_number(total_AM_primitives)=5
    
      CALL read_modifiers(line,length,pos)
    
    else if (ch.EQ.'6') then
  
      write(*,*)'Read Moire primitive'
    
      AM_primitive_number(total_AM_primitives)=6
    
      CALL read_modifiers(line,length,pos)
    
    else if (ch.EQ.'7') then
  
      write(*,*)'Read thermal primitive'
    
      AM_primitive_number(total_AM_primitives)=7
    
      CALL read_modifiers(line,length,pos)
    
    end if
    
  end if

end do

if (line(length:length).NE.'%') then
  write(*,*)'Error in Aperture macro:',trim(line)
  write(*,*)'ch:',line(length:length)
  STOP
end if

RETURN

END SUBROUTINE read_Aperture_Macro
!
! ____________________________________________________________
!
!
SUBROUTINE read_modifiers(line,length,pos)

! read a list of modifiers, sperated by commas until a '*' (end of block) is reached

USE gerber
IMPLICIT NONE

character*256 :: line
integer   :: length,pos

! local variables

integer :: modifier
integer :: i
character :: ch

! START

modifier=0
i=0

do while (pos.LT.length)

  pos=pos+1
  ch=line(pos:pos)
  
  if (ch.EQ.',') then
  
    modifier=modifier+1
    i=0
    
  else if (ch.EQ.'*') then
  
    exit
    
  else
  
    i=i+1
! replace 'x' and 'X' with *
    if ( (ch.EQ.'x').OR.(ch.EQ.'X') ) ch='*'
    AM_modifiers(total_AM_primitives,modifier)(i:i)=ch
    
  end if
  
end do

n_AM_modifiers(total_AM_primitives)=modifier

write(*,*)'Number of modifiers read:',modifier

do i=1,modifier
  write(*,*)i,'  :',trim( AM_modifiers(total_AM_primitives,i) )
end do

RETURN

END SUBROUTINE read_modifiers
!
! ____________________________________________________________
!
!
SUBROUTINE read_variable_definition(line,length,pos)

! read the variable definition until a '*' (end of block) is reached

USE gerber
IMPLICIT NONE

character*256 :: line
integer   :: length,pos

! local variables

integer :: modifier
integer :: i
character :: ch

! START

modifier=0
i=0
! check for '='

ch=line(pos:pos)
if (ch.NE.'=') then
  write(*,*)"Error reading variable definition. Expecting '=', found:",ch,' at position',pos
  STOP
end if

! read a string and assing to modifier 1
modifier=1

do while (pos.LT.length)

  pos=pos+1
  ch=line(pos:pos)
  
  if (ch.EQ.'*') then
  
    exit
    
  else
  
    i=i+1
! replace 'x' and 'X' with *
    if ( (ch.EQ.'x').OR.(ch.EQ.'X') ) ch='*'
    AM_modifiers(total_AM_primitives,modifier)(i:i)=ch
    
  end if
  
end do

n_AM_modifiers(total_AM_primitives)=modifier

write(*,*)'Variable definition string:',trim( AM_modifiers(total_AM_primitives,modifier) )

RETURN

END SUBROUTINE read_variable_definition
