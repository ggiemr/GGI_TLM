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
! Name rectangular_cavity
!     
!
! Description
!     
!
! Comments:
!      
!
! History
!
!     started dd/mm/yy CJS
!

PROGRAM rectangular_cavity

USE constants

IMPLICIT NONE

! local_variables

  real*8 lx,ly,lz,epsr,mur,c
  integer m_max,n_max,p_max
  integer m_min,n_min,p_min
  integer m,n,p
  
  integer nf,i,ii
  real*8  fswap
  integer iswap(3)
  
  real*8,allocatable :: f(:)
  integer, allocatable :: mnp(:,:)
  

! function_types

! START

  write(*,*)'CALLED: rectangular_cavity'

  write(*,*)'Enter cavity dimensions in m'
  read(*,*)lx,ly,lz
  write(*,*)'Enter relative permittivity'
  read(*,*)epsr
  write(*,*)'Enter relative permeability'
  read(*,*)mur
  write(*,*)'Enter minimum m, n and p'
  read(*,*)m_min,n_min,p_min
  write(*,*)'Enter maxmimum m, n and p'
  read(*,*)m_max,n_max,p_max

  nf=(n_max-n_min+1)*(m_max-m_min+1)*(p_max-p_min+1)

  allocate(f(1:nf))
  allocate(mnp(1:nf,1:3))
  
  c=c0/sqrt(epsr*mur)

  i=0

  do m=m_min,m_max
    do n=n_min,n_max
      do p=p_min,p_max
      
        i=i+1
	mnp(i,1)=m
	mnp(i,2)=n
	mnp(i,3)=p
	if ((m+n.eq.0).OR.(m+p.eq.0).OR.(n+p.eq.0)) then
! no mode if more than one of m,n or p is zero
	  f(i)=0d0
	else
	  f(i)=(c/2d0)*sqrt((m/lx)**2+(n/ly)**2+(p/lz)**2)
	end if
      end do
    end do
  end do

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

  open(unit=10,file='rectangular_cavity_modes')

  do i=1,nf
    if (f(i).ne.0d0) then
      write(10,1000)mnp(i,1:3),f(i)
    end if
  end do
1000 format(3I5,E16.8)
  
  close(unit=10)

  deallocate(f)
  deallocate(mnp)

  write(*,*)'FINISHED: rectangular_cavity'
  
END 

  
  
  
  
