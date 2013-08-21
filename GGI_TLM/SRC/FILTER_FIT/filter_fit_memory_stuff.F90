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
! SUBROUTINE deallocate_filter_fit_memory
!
! NAME
!     deallocate_filter_fit_memory
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/12/2012 CJS
!
!
SUBROUTINE deallocate_filter_fit_memory

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_testing_data
USE FF_filters
USE FF_file_stuff

IMPLICIT NONE
 
! local variables

  integer i

! START

! deallocate input data stuff
  if (allocated( frequency )) DEALLOCATE ( frequency )
  if (allocated( normalised_frequency )) DEALLOCATE ( normalised_frequency )
  if (allocated( s )) DEALLOCATE ( s )
  if (allocated( value ))    DEALLOCATE ( value )

! deallocate testing stuff  
  if (allocated( testing_frequency) ) DEALLOCATE( testing_frequency )
  if (allocated( testing_s) ) DEALLOCATE( testing_s )
  
! deallocates filters
  if ( allocated( filter_S ) ) then
    
    do i=1,n_functions
      CALL deallocate_Sfilter( filter_S(i) )
    end do
    
    DEALLOCATE( filter_S )
  
  end if
  
  if ( allocated( filter_sigma ) ) DEALLOCATE( filter_sigma )


END SUBROUTINE deallocate_filter_fit_memory
