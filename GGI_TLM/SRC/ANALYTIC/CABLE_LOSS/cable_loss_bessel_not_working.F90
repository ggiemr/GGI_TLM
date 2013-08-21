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

real*8 fmin,fmax,fstep,f
integer floop,nf

real*8 delta,q

complex*16 qc

complex*16 J0,J0p

real*8 ber,bei,berp,beip

real*8 rdc,lidc

real*8 r,li

! function variables

real*8 derivative_besselJ

! START

  write(*,*)'Calculate the frequency dependent impedance due to skin depth of a cylindrical wire'

  write(*,*)'Enter the wire radius (m)'
  read(*,*)rw

  write(*,*)'Enter the wire conductivity (S/m)'
  read(*,*)sigma
    
  write(*,*)'Enter frequency minimum, maximum '
  read(*,*)fmin,fmax
    
  write(*,*)'Enter the number of frequency steps '
  read(*,*)nf
  
  fstep=(fmax-fmin)/(nf-1)
  
  open(unit=10,file='cable_loss.fout')
    
  rdc=1d0/(sigma*pi*rw*rw)
    
  lidc=mu0/(8d0*pi)
    
! loop over frequency
  do floop=1,nf
  
    f=fmin+(floop-1)*fstep

    delta=1d0/sqrt(pi*f*mu0*sigma)
    
    q=sqrt(2d0)*rw/delta
    qc=sqrt(-2d0*j)*rw/delta
    
    J0 =besjn (0,qc)
    J0p=derivative_besselJ(0,qc)
    
    ber =dble(J0)
    bei =imag(J0)
    berp=dble(J0p)
    beip=imag(J0p)
    
    r=rdc*(q/2d0)*(ber*beip-bei*berp)/(beip*beip+berp*berp)    
    
    li=lidc*(4d0/q)*(bei*beip+ber*berp)/(beip*beip+berp*berp)    
            
    write(10,8000)f,r,li

8000 format(3E16.6)
  
  end do ! next frequency
  
  close(unit=10)
  
END PROGRAM cable_loss


