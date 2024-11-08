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
! Tried to make this work for arbitrary angles and polarisation. 
! NEEDS THOROUGH TESTING CJS 27/1/2015
!
!
! HISTORY
!
!     started 12/01/10 CJS
!     Attempt to get it working for arbitrary incicent angles CJS 27/1/2015
!
  PROGRAM Mlayer

USE constants
USE filter_types
USE filter_operators
USE filter_functions

USE Mlayer_file_module
USE Mlayer_module

IMPLICIT NONE

! local variables

  integer n_frequencies,frequency_loop

! START

  write(*,*)'Calculate reflection and transmission of multiple material layers '
  write(*,*)'and impedance boundariesas a function of frequency,'
  write(*,*)'angle of incidence and polarisation'
  write(*,*)'Note:'
  write(*,*)'The material layers are normal to the z direction'
  write(*,*)'Wave propagation is in the xz plane'
  write(*,*)'There is no field variation in the y direction (normal to the plane of propagation)'
  write(*,*)'Angle of incidence is measured from the normal to the layers'
  write(*,*)'i.e. the angle between the k vector and the z direction'
  write(*,*)'TE(z) polarisation: Ex, Hy, Ez'
  write(*,*)'TM(z) polarisation: Hx, Ey, Hz'
  
  call open_files()
  
  call read_Mlayer_input_file()
  
  overflow_error=.FALSE.
  
  if (logf_flag) then
    n_frequencies=nf
    logfmin=log10(fmin)
    logfmax=log10(fmax)
    logfstep=(logfmax-logfmin)/dble(nf-1)
  else
    n_frequencies=int( (fmax-fmin)/fstep )+1
  end if
  
  do frequency_loop=1,n_frequencies

    if (logf_flag) then
      logf=logfmin+(frequency_loop-1)*logfstep
      f=10d0**logf
    else
      f=fmin+(frequency_loop-1)*fstep    
    end if
    
    w=2d0*pi*f
  
    call calc_Sparams()
  
  end do
  
  call deallocate_memory()
  
  call close_files()
  
  if (overflow_error) then
    write(*,*)'WARNING: Overflow error occurred in Mlayer'
  end if

  END PROGRAM Mlayer
