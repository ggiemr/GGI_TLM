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
!  SUBROUTINE allocate_memory
!  SUBROUTINE deallocate_memory


! NAME
!     SUBROUTINE  allocate_memory
!
! DESCRIPTION
! 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
  SUBROUTINE allocate_memory()

USE constants
USE filter_types
USE filter_operators
USE filter_functions

USE Mlayer_file_module
USE Mlayer_module

! START

  write(info_file_unit,*)'Allocating memory'
  
  allocate(layer_type(1:m))
  allocate(layer_thickness(1:m))
  allocate(material_list(1:m))
  allocate(thin_layer_list(1:m))
  allocate(ABCD(1:2,1:2,1:m))

  return

  END SUBROUTINE allocate_memory
  
! NAME
!     SUBROUTINE  deallocate_memory
!
! DESCRIPTION
! 
!	
!     
! COMMENTS
!     
! 
!
!
! HISTORY
!
!     started 12/01/10 CJS
!
!
  SUBROUTINE deallocate_memory()

USE constants
USE filter_types
USE filter_operators
USE filter_functions

USE Mlayer_file_module
USE Mlayer_module

! START

  write(info_file_unit,*)'Deallocating memory'

  deallocate(layer_type)
  deallocate(layer_thickness)
  deallocate(material_list)
  deallocate(thin_layer_list)
  deallocate(ABCD)


  return

  END SUBROUTINE deallocate_memory
  
