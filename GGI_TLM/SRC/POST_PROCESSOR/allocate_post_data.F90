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
! SUBROUTINE allocate_post_data
! SUBROUTINE deallocate_post_data
!
! NAME
!    allocate_post_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 18/1/2013 CJS
!
!
SUBROUTINE allocate_post_data

USE post_process

IMPLICIT NONE

! local variables

! START

  if (n_functions_of_time.NE.0) then
    ALLOCATE( function_of_time(1:n_functions_of_time) )
  end if
  
  if (n_functions_of_frequency.NE.0) then
    ALLOCATE( function_of_frequency(1:n_functions_of_frequency) )
  end if
  
  RETURN
  
  
END SUBROUTINE allocate_post_data
!
! NAME
!    deallocate_post_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 18/1/2013 CJS
!
!
SUBROUTINE deallocate_post_data

USE post_process

IMPLICIT NONE

! local variables

  integer :: function

! START

  if (allocated( function_of_time ) ) then
  
    do function=1,n_functions_of_time
    
      if (allocated( function_of_time(function)%time )) DEALLOCATE( function_of_time(function)%time )
      if (allocated( function_of_time(function)%value )) DEALLOCATE( function_of_time(function)%value )
      
    end do ! next function
    
    DEALLOCATE( function_of_time )
    
  end if
  
  if (allocated( function_of_frequency ) ) then
  
    do function=1,n_functions_of_frequency
    
      if (allocated( function_of_frequency(function)%frequency )) DEALLOCATE( function_of_frequency(function)%frequency )
      if (allocated( function_of_frequency(function)%value )) DEALLOCATE( function_of_frequency(function)%value )
      if (allocated( function_of_frequency(function)%magnitude )) DEALLOCATE( function_of_frequency(function)%magnitude )
      if (allocated( function_of_frequency(function)%phase )) DEALLOCATE( function_of_frequency(function)%phase )
      if (allocated( function_of_frequency(function)%dB )) DEALLOCATE( function_of_frequency(function)%dB )
      
    end do ! next function
    
    DEALLOCATE( function_of_frequency )
    
  end if
  
  
  RETURN
  
  
END SUBROUTINE deallocate_post_data
