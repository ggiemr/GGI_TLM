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
!	convert_to_lower_case
!
! DESCRIPTION
!     returns the input string with all characters converted to lower case
!
! HISTORY
!
!     started 10/02/09 CJS
!
! COMMENTS
!     
subroutine convert_to_lower_case(input_string,input_string_length)

! variables passed to subroutine

  integer input_string_length
  character*(input_string_length) input_string
  
! local variables

  integer i
  character*(input_string_length) output_string
  integer na,nCA,nz,nCZ,diff
  integer c
  
! START

  output_string=input_string

  na=IACHAR('a')
  nCA=IACHAR('A')
  nz=IACHAR('z')
  nCZ=IACHAR('Z')
  
  diff=na-nCA
  
  do i=1,input_string_length
  
    c=IACHAR(input_string(i:i))
    if ((c.ge.nCA).AND.(c.le.nCZ)) then
      c=c+diff
      output_string(i:i)=ACHAR(c)
    end if
    
  end do
  
  input_string=output_string
  
  return
  
end subroutine convert_to_lower_case
!
! NAME
!	convert_to_upper_case
!
! DESCRIPTION
!     returns the input string with all characters converted to upper case
!
! HISTORY
!
!     started 26/02/09 CJS
!
! COMMENTS
!     
subroutine convert_to_upper_case(input_string,input_string_length)

! variables passed to subroutine

  integer input_string_length
  character*(input_string_length) input_string
  
! local variables

  integer i
  character*(input_string_length) output_string
  integer na,nCA,nz,nCZ,diff
  integer c
  
! START

  output_string=input_string

  na=IACHAR('a')
  nCA=IACHAR('A')
  nz=IACHAR('z')
  nCZ=IACHAR('Z')
  
  diff=nCA-na
  
  do i=1,input_string_length
  
    c=IACHAR(input_string(i:i))
    if ((c.ge.na).AND.(c.le.nz)) then
      c=c+diff
      output_string(i:i)=ACHAR(c)
    end if
    
  end do
  
  input_string=output_string
  
  return
  
end subroutine convert_to_upper_case
