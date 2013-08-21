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
  PROGRAM RCS_cylinder
    
! calculate the RCS of a finite (or SW of infinite) circular cylinder

USE constants

IMPLICIT NONE

integer max_terms
parameter (max_terms=50)

real*8 error
parameter (error=1D-6)

real*8 a,l

real*8 fmin,fmax,fstep,f,lambda,w
integer floop,nf

real*8 tmin,tmax,tstep,t,t_rad
integer tloop,nt

integer n,n_terms,epsilon

real*8 beta,beta_a

complex*16 term_TE(0:max_terms)
complex*16 term_TM(0:max_terms)

real*8 sum_TE,sum_TM
complex*16 Csum_TE,Csum_TM
real*8 sigma_TE,sigma_TM

logical convert_3D

! function variables

real*8 derivative_besselJ
real*8 derivative_besselY

! START

  write(*,*)'Calculation of the RCS of circular cylinders'
  write(*,*)'Enter the cylinder radius (m)'
  read(*,*)a 
  
  write(*,*)'Enter the cylinder (m) length or 0 for calculation of SW of infinite cylinder'
  read(*,*)l
  
  write(*,*)'Enter frequency minimum, maximum and step'
  read(*,*)fmin,fmax,fstep
  
  write(*,*)'Enter angle minimum, maximum and step'
  read(*,*)tmin,tmax,tstep
  
  convert_3D=.FALSE.
  if (l.ne.0d0) convert_3D=.TRUE.
  
  nf=1+NINT((fmax-fmin)/fstep)
  
  nt=1+NINT((tmax-tmin)/tstep)
  
  open(unit=10,file='RCS_cylinder.dat')
  
  if (convert_3D) then
    write(*,*)'Calculating Radar Cross Section (dBsm)'
    write(10,*)'# frequency(Hz)  angle(degrees) sigma_3D_TE(dBsm)  sigma_3D_TM(dBsm)  n_terms'
    write(10,*)'#                               (E perpendicular)   (E Parallel)   '
  else
    write(*,*)'Calculating Scattering Width (dBm)'
    write(10,*)'# frequency(Hz)  angle(degrees) sigma_2D_TE(dBm)  sigma_2D_TM(dBm)  n_terms'
    write(10,*)'#                               (E perpendicular)   (E Parallel)   '
  end if
  
! loop over frequency
  do floop=1,nf
  
    f=fmin+(floop-1)*fstep
    lambda=c0/f
    w=2d0*pi*f
    beta=w/c0
    beta_a=beta*a
    
! calculaTM the series of terms of the expansion
    sum_TM=0d0
    sum_TE=0d0
    
    do n=0,max_terms
      
      if (n.eq.0) then
        epsilon=1
      else
        epsilon=2
      end if
      
      term_TM(n)=epsilon*besjn (n,beta_a)/	&
                 (besjn (n,beta_a)-j*besyn (n,beta_a))
      sum_TM=sum_TM+(abs(term_TM(n)))**2
      
      term_TE(n)=epsilon*derivative_besselJ (n,beta_a)/	&
                 (derivative_besselJ (n,beta_a)-j*derivative_besselY (n,beta_a))
      sum_TE=sum_TE+(abs(term_TE(n)))**2
      
!      write(*,*)n,real(term_TM(n)),imag(term_TM(n)),real(term_TE(n)),imag(term_TE(n))
      
      if ( ((abs(term_TM(n)))**2.lt.error*sum_TM) .AND.	&
           ((abs(term_TE(n)))**2.lt.error*sum_TE) ) then
! assume that the series has converged
        n_terms=n
	GOTO 1000
      end if
    
    end do
    
    write(*,*)'Warning, max_terms reached in expansion, f=',f
    n_terms=max_terms

1000 CONTINUE
    
! loop over angle
    do tloop=1,nt
    
      t=(tmin+(tloop-1)*tstep)   
      t_rad=t*pi/180d0   ! note convert to radians here...

! assemble the series solution for this angle      
      Csum_TM=(0d0,0d0)
      Csum_TE=(0d0,0d0)
      do n=0,n_terms
	
        Csum_TM=Csum_TM+term_TM(n)*cos(n*t_rad) 
      
        Csum_TE=Csum_TE+term_TE(n)*cos(n*t_rad) 
  
      end do ! next term in the series
      
      sigma_TM=((abs(Csum_TM))**2)*2.0*lambda/pi
      sigma_TE=((abs(Csum_TE))**2)*2.0*lambda/pi
       
      if (convert_3D) then ! convert from scatTMring width to RCS using the cylinder length
        sigma_TM=sigma_TM*2d0*l*l/lambda
        sigma_TE=sigma_TE*2d0*l*l/lambda
      end if
            
      write(10,8000)f,t,10d0*log10(sigma_TE),10d0*log10(sigma_TM),n_terms
!      write(10,8000)f,t,sigma_TE/lambda,sigma_TM/lambda,n_terms
8000  format(4E16.6,I10)
      
    end do ! next angle
  
  end do ! next frequency
  
  close(unit=10)
  
END PROGRAM RCS_cylinder


