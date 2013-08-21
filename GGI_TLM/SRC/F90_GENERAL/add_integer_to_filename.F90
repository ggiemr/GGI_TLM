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
!SUBROUTINE add_integer_to_filename
!
! NAME
!     SUBROUTINE add_integer_to_filename
!
! DESCRIPTION
!     
!  add an integer to a given filename
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
  SUBROUTINE add_integer_to_filename(filename,order,opfilename)
!
    character*(*) filename,opfilename
    integer order
    character chorder
    character*20 stringorder
    integer namlen

    integer n_char_order,i,string_pos
    integer digit,number,n_digits
    integer mask

! START

    filename_len=len(filename)
       
    if (order.eq.0) then
      n_digits=1
      stringorder(1:1)='0'
    else   
       
      number=order
      n_digits=int(log10(dble(number)))+1
    
      mask=1
    
      do i=1,n_digits

        digit=mod(number,mask*10)/mask
        string_pos=n_digits-i+1      
        stringorder(string_pos:string_pos)=char(digit+ichar('0'))
        number=number-mask*digit
        mask=mask*10
      
      end do ! next digit
      
    end if
    
    opfilename=trim(filename)//"."//stringorder(1:n_digits)
    
    RETURN
    
    END SUBROUTINE add_integer_to_filename
