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
  PROGRAM cable_loss
    
! calculate the frequency dependent impedance due to skin depth of a cylindrical wire

! See C. Paul, "Analysis of multiconductor transmission lines" Wiley, pp164-169

USE constants

IMPLICIT NONE

  real*8 rw,sigma

  real*8        :: f,fmin,fmax,fstep
  real*8        :: log_fmin,log_fmax,log_fstep,log_f
  integer       :: n_frequencies
  integer       :: frequency_loop
  
  character*3	:: freq_range_type

  real*8 delta,f0

  complex*16 Zc

  real*8 rdc

  real*8 r,li

! function variables

real*8 derivative_besselJ

! START

  write(*,*)'Calculate the frequency dependent impedance due to skin depth of a cylindrical wire'

  write(*,*)'Enter the wire radius (m)'
  read(*,*)rw

  write(*,*)'Enter the wire conductivity (S/m)'
  read(*,*)sigma
    
  write(*,*)'Enter the frequency range type (log or lin)'
  read(*,'(A3)')freq_range_type
  
  if ( (freq_range_type.NE.'log').AND.(freq_range_type.NE.'lin') ) then
    write(*,*)"Frequency range type should be 'log' or 'lin'"
    STOP
  end if
  
  write(*,*)'Enter minimum frequency, fmin'
  read(*,*)fmin
  write(*,*)'Enter maximum frequency, fmax'
  read(*,*)fmax
  write(*,*)'Enter the number of frequencies'
  read(*,*)n_frequencies
  
  if (freq_range_type.EQ.'log') then
  
    log_fmin=log10(fmin)
    log_fmax=log10(fmax)
    log_fstep=(log_fmax-log_fmin)/dble(n_frequencies-1)
  
  else if (freq_range_type.EQ.'lin') then
  
    fstep=(fmax-fmin)/dble(n_frequencies-1)
  
  end if
  
  open(unit=10,file='cable_loss.fout')
    
  rdc=1d0/(sigma*pi*rw*rw)
  
  f0=4d0/(rw*rw*pi*mu0*sigma)
  
  write(*,*)'f0=',f0
    
! loop over frequency
  
  do frequency_loop=1,n_frequencies
  
    if (freq_range_type.EQ.'log') then
    
      log_f=log_fmin+(frequency_loop-1)*log_fstep

      f=10D0**log_f
  
    else if (freq_range_type.EQ.'lin') then
    
      f=fmin+(frequency_loop-1)*fstep
  
    end if
    
!    delta=1d0/sqrt(pi*f*mu0*sigma)  ! not used

    if (f.lt.f0) then
      Zc=rdc*(1d0+j*f/f0)
    else
      Zc=rdc*sqrt(f/f0)*(1d0+j)
    end if
    
    r=dble(Zc)
    
    li=dimag(Zc)    
            
    write(10,8000)f,r,li

8000 format(3E16.6)
  
  end do ! next frequency
  
  close(unit=10)
  
END PROGRAM cable_loss


