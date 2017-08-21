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
! Name rectangular_waveguide_modes
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

PROGRAM rectangular_waveguide_modes

IMPLICIT NONE

! local_variables

  real*8 a,b,epsr,mur,c
  integer m_max,n_max
  integer m_min,n_min
  integer m,n
  
  integer nf,i,ii
  real*8  fswap
  integer iswap(2)
  
  real*8,allocatable :: f(:)
  integer, allocatable :: mn(:,:)
  
  real*8,parameter :: pi=3.1415926535d0
  real*8,parameter :: c0=2.998D8
  
! function_types

! START

  write(*,*)'CALLED: rectangular_waveguide_modes'

  write(*,*)'Enter waveguide width, a, and height, b, in metres'
  read(*,*)a,b
  write(*,*)'Enter relative permittivity'
  read(*,*)epsr
  write(*,*)'Enter relative permeability'
  read(*,*)mur
  write(*,*)'Enter minimum m and n'
  read(*,*)m_min,n_min
  write(*,*)'Enter maxmimum m and n '
  read(*,*)m_max,n_max

  nf=(n_max-n_min+1)*(m_max-m_min+1)

  allocate(f(1:nf))
  allocate(mn(1:nf,1:2))
  
  c=c0/sqrt(epsr*mur)

  i=0

  do m=m_min,m_max
    do n=n_min,n_max
    
      write(*,*)'Cutoff of mode m=',m,' n=',n
      
      i=i+1
      mn(i,1)=m
      mn(i,2)=n
      if (m+n.eq.0) then
! no mode if more than one of m or n  is zero
	f(i)=0d0
        write(*,*)'No mode exists'
      else
	f(i)=(c/(2d0*pi))*sqrt((m*pi/a)**2+(n*pi/b)**2)
        write(*,'(A2,ES16.6)')'f=',f(i)
      end if
        
    end do
  end do

! sort mode frequencies (not a very efficient way of sorting here...)

  do i=1,nf-1
! find the ith smallest in this loop
    do ii=i+1,nf
    
      if (f(ii).lt.f(i)) then
! swap
        fswap=f(i)
	iswap(1:2)=mn(i,1:2)
	f(i)=f(ii)
	mn(i,1:2)=mn(ii,1:2)
	f(ii)=fswap
	mn(ii,1:2)=iswap(1:2)
        
      end if
    
    end do
  end do

  open(unit=10,file='rectangular_waveguide_modes_modes')
  
  write(*,*)'Modes ordered according to cutoff frequency'

  do i=1,nf
    if (f(i).ne.0d0) then
      write(10,1000)mn(i,1:2),f(i)
      write(*,1000)mn(i,1:2),f(i)
    end if
  end do
1000 format(2I5,ES16.8)
  
  close(unit=10)

  deallocate(f)
  deallocate(mn)

  write(*,*)'FINISHED: rectangular_waveguide_modes'
  
END 

  
