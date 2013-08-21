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
!     PROGRAM Mlayer
!
! DESCRIPTION
! calculate the S parameters of multiple layers comprising 
! bulk materials and impedance boundary layers
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
  PROGRAM Mlayer

USE constants
USE filter_types
USE filter_operators
USE filter_functions

USE Mlayer_file_module
USE Mlayer_module

! local variables

  integer n_frequencies,frequency_loop

! START

  write(*,*)'Calculate reflection and transmission of layers in the yz plane'
  write(*,*)'as a function of frequency, angle of incidence'
  write(*,*)'and polarisation'
  write(*,*)'Note:'
  write(*,*)'1. Angle of incidence is measured from normal'
  write(*,*)'TE(z) polarisation: Ex, Ey, Hz'
  write(*,*)'TM(z) polarisation: Ez, Hx, Hy'
  
  call open_files()
  
  call read_Mlayer_input_file()
  
  overflow_error=.FALSE.
  
  n_frequencies=int( (fmax-fmin)/fstep )+1

  do frequency_loop=1,n_frequencies

    f=fmin+(frequency_loop-1)*fstep
  
    w=2d0*pi*f
  
    call calc_Sparams()
  
  end do
  
  call deallocate_memory()
  
  call close_files()
  
  if (overflow_error) then
    write(*,*)'WARNING: Overflow error occurred in Mlayer'
  end if

  END PROGRAM Mlayer
