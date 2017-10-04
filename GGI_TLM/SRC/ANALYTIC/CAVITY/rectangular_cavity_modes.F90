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
! Name rectangular_cavity_modes
!     
!
! Description
!     Calculate the resonant frequencies of a rectangular cavity
!
! Comments:
!      
!
! History
!
!     started 4/10/2017 CJS
!

PROGRAM rectangular_cavity_modes

IMPLICIT NONE

! local_variables

  real*8  w,d,h
  real*8  epsr,mur,c
  integer m_max,n_max,p_max
  integer m_min,n_min,p_min
  integer m,n,p
  integer nzeroindex
  
  integer nf,i,ii
  real*8  fswap
  integer iswap(2)
  
  real*8,allocatable :: f(:)
  integer, allocatable :: mnp(:,:)
  
  real*8,parameter :: pi=3.1415926535d0
  real*8,parameter :: c0=2.998D8
  
! function_types

! START

  write(*,*)'CALLED: rectangular_cavity_modes'

  write(*,*)'Enter waveguide width (x), w, depth (y), d, and height (z), h, in metres'
  read(*,*)w,h,d
  write(*,*)'Enter relative permittivity'
  read(*,*)epsr
  write(*,*)'Enter relative permeability'
  read(*,*)mur
  write(*,*)'Enter minimum mode index in each direction, m, n and p'
  read(*,*)m_min,n_min,p_min
  write(*,*)'Enter maxmimum mode index in each direction, m, n and p'
  read(*,*)m_max,n_max,p_max

  nf=(n_max-n_min+1)*(m_max-m_min+1)*(p_max-p_min+1)

  allocate(f(1:nf))
  allocate(mnp(1:nf,1:3))
  
  c=c0/sqrt(epsr*mur)

  i=0

  do m=m_min,m_max
    do n=n_min,n_max
      do p=p_min,p_max
    
        write(*,*)'Cutoff of mode m=',m,' n=',n,' p=',p
      
        i=i+1
        mnp(i,1)=m
        mnp(i,2)=n
        mnp(i,3)=p
        
! if more than one index m,n,p is zero then no mode exists
        nzeroindex=0
        if (m.EQ.0) nzeroindex=nzeroindex+1
        if (n.EQ.0) nzeroindex=nzeroindex+1
        if (p.EQ.0) nzeroindex=nzeroindex+1
        
        if (nzeroindex.GT.1) then
        
! no mode if more than one of m or n  is zero
	  f(i)=0d0
          write(*,*)'No mode exists'
          
        else
        
	  f(i)=(c/(2d0*pi))*sqrt((m*pi/w)**2+(n*pi/d)**2+(p*pi/h)**2)
          write(*,'(A2,ES16.6)')'f=',f(i)
          
        end if
        
      end do   ! next p
    end do   ! next n
  end do   !next m

! sort mode frequencies (not a very efficient way of sorting here...)

  do i=1,nf-1
! find the ith smallest in this loop
    do ii=i+1,nf
    
      if (f(ii).lt.f(i)) then
! swap
        fswap=f(i)
	iswap(1:3)=mnp(i,1:3)
	f(i)=f(ii)
	mnp(i,1:3)=mnp(ii,1:3)
	f(ii)=fswap
	mnp(ii,1:3)=iswap(1:3)
        
      end if
    
    end do
  end do

  open(unit=10,file='rectangular_cavity_modes_modes')
  
  write(*,*)'Modes ordered according to cutoff frequency'

  do i=1,nf
    if (f(i).ne.0d0) then
      write(10,1000)mnp(i,1:3),f(i)
      write(*,1000)mnp(i,1:3),f(i)
    end if
  end do
1000 format(3I5,ES16.8)
  
  close(unit=10)

  deallocate(f)
  deallocate(mnp)

  write(*,*)'FINISHED: rectangular_cavity_modes'
  
END 

  
