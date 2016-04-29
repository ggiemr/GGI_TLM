! SUBROUTINE wigner_function
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
!    wigner_function
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 23/10/2014 CJS
!     attempt to improve the numerics of the Wigner function propagation
!     attempt to improve the speed of the process 20/1/2016 CJS: reorder arrays, pre-calculate exponential functions, single precision
!
SUBROUTINE wigner_function

USE post_process
USE file_information
USE Constants

IMPLICIT NONE

! local variables

integer			:: nx,ny
integer			:: np,npx,npy

  integer		:: max_col
  integer		:: n_lines
  real*4		:: value(4)
  integer		:: step(0:4)
  integer		:: n(0:4)
  real*4		:: data_line(1:6)
  integer		:: loop
  integer		:: num,den
  
  real*4,allocatable    :: x_values(:)
  real*4,allocatable    :: y_values(:)
  
  real*4,allocatable    :: sx_values(:)
  real*4,allocatable    :: sy_values(:)
  

complex,allocatable	:: Cov(:,:,:,:)
complex,allocatable	:: WF(:,:,:,:)

complex,allocatable	:: WF2(:,:,:,:)
complex,allocatable	:: Cov2(:,:,:,:)
integer,allocatable	:: Nsamples_WF2(:,:,:,:)

complex,allocatable	:: exp_function(:,:,:,:)

real*4			:: frequency
real*4			:: dx,dy,dA

character*256		:: filename

real*4			:: k,w

real*4			:: x1,xmin,xmax,x2
real*4			:: y1,ymin,ymax,y2
real*4			:: sx2,sx,sy2,sy
real*4			:: px,py,dp,kdpA
real*4			:: pmin,pmax,mag_p

real*4			:: z

complex		:: sum

real*4			:: re,im,x0
integer			:: element_x1,element_y1,element_x2,element_y2

integer			:: ix1,ix2,iy1,iy2,isx,isy,ipx,ipy
integer			:: array_sx,array_sy

integer			:: i,ii,i1,i2,i3,i4

character :: ch
logical	:: inc_evanescent


! START
  
! READ COVARIANCE MATRIX (correlation matrix)

! NEED TO WORK OUT WHETHER WE HAVE 1D OR 2D DATA AND THE DIMENSIONS OF THE DATA SET (NX1,NY1,NX2,NY2)

! ENTER OPTIONS FOR WIGNER FUNCTION CALCULATION

  write(*,*)'Enter the filename for the input covariance data'
  read(*,'(A)')filename
  open(unit=local_file_unit,file=filename)
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  max_col=6
   
! work out the size of each dimension in the file

  n_lines=0
  step(:)=0
  value(:)=0.0
   
10 CONTINUE

    read(local_file_unit,*,end=100)(data_line(i),i=1,max_col)
    
    n_lines=n_lines+1
    
    if (n_lines.EQ.1) then
    
! read origin point for multi-column data
      value(1:4)=data_line(1:4)
      
    else

! work out how the x and y values change in the file
      do i=1,4
        if ( (data_line(i).NE.value(i)).AND.(step(i).eq.0) ) step(i)=n_lines-1
      end do
      
    end if
    
    GOTO 10
    
100 CONTINUE
  
  step(0)=n_lines
  
! work out the number of elements in each dimensions
  
  do i=0,4
    write(*,'(A,I2,A,I10)')'Step ',i,'=',step(i)
  end do
  
  n(:)=1 ! default dimension is 1
  
  do i=1,4
  
    if (step(i).NE.0) then
! there is some data in this dimension
      den=step(i)
! find the next largest step value
      num=999999999
      do loop=0,4
        if ( (step(loop).GT.den).AND.(step(loop).LT.num) ) then
          num=step(loop)
        end if
      end do
      n(i)=num/den
      
  else
! there is no data in this dimension

    n(i)=1
    
  end if
  
  end do
  
  do i=1,4
    write(*,'(A,I2,A,I10)')'n ',i,'=',n(i)
  end do
  
  write(*,*)'nx1=',n(1)
  write(*,*)'ny1=',n(2)
  write(*,*)'nx2=',n(3)
  write(*,*)'ny2=',n(4)
  
! check
  if ( (n(1).NE.n(3)).OR.(n(2).NE.n(4)) ) then
    write(*,*)'Error reading correlation matrix in Wigner function process'
    write(*,*)' (nx1.NE.nx2).OR.(ny1.NE.ny2) '
    STOP
  end if
      
! get the plotting dimensions
  nx=n(1)
  ny=n(2)
  
  write(*,*)'Enter the number of samples in p required'
  read(*,*)np
  write(record_user_inputs_unit,*)np,'  ! number of samples in p required'
  
  npx=np
  npy=np
  
! limit the p values if the correlation results from 1D field scan data  
  if (n(1).EQ.1)  npx=1
  if (n(2).EQ.1)  npy=1
  

  write(*,*)'ALLOCATE MEMORY AND READ COVARIANCE DATA FILE'
  
  ALLOCATE(Cov(1:nx,1:nx,1:ny,1:ny))
  
  ALLOCATE(WF(1:np,1:np,1:2*nx-1,1:2*ny-1))
  
  ALLOCATE(Cov2(1:nx,1:nx,1:ny,1:ny))
  
  ALLOCATE(WF2(1:np,1:np,1:2*nx-1,1:2*ny-1))
  ALLOCATE(Nsamples_WF2(1:np,1:np,1:2*nx-1,1:2*ny-1))
  
  ALLOCATE( x_values(1:nx) )
  ALLOCATE( y_values(1:ny) )
  
  ALLOCATE( sx_values(1:2*nx+1) )
  ALLOCATE( sy_values(1:2*ny+1) )

  rewind(local_file_unit)
  
  x_values(1:nx)=0.0
  y_values(1:ny)=0.0
  
  do ix1=1,nx
    do iy1=1,ny
      do ix2=1,nx
        do iy2=1,ny
    
          read(local_file_unit,*)x1,y1,x2,y2,re,im

! We could check the file format here maybe???
! note 21/1/2016: reordered array
          x_values(ix1)=x1
          y_values(iy1)=y1
          Cov(ix1,ix2,iy1,iy2)=cmplx(re,im)
      
        end do     
      end do     
    end do     
  end do
  
  close(unit=local_file_unit)
  
  dx=x_values(2)-x_values(1)
  dy=y_values(2)-y_values(1)
  
  xmax=nx*dx/2.0
  xmin=-xmax
  do ix1=1,2*nx-1
    sx_values(ix1)=xmin+ix1*dx/2.0
  end do
  
  ymax=ny*dy/2.0
  ymin=-ymax
  do iy1=1,2*ny-1
    sy_values(iy1)=ymin+iy1*dy/2.0
  end do


! we now read x and y values from the correlation input file  
!  write(*,*)'Enter the spatial resolution of the samples'
!  read(*,*)dl
!  write(record_user_inputs_unit,*)dl,'  ! dl'
  
  write(*,*)'Enter the frequency'
  read(*,*)frequency
  write(record_user_inputs_unit,*)frequency,'  ! frequency'
  
  write(*,*)'Enter pmax (note range of p is -pmax to pmax, pmax<=1 excludes evanescent waves)'
  read(*,*)pmax
  write(record_user_inputs_unit,*)pmax,'  !  pmax  (note range of p is -pmax to pmax, pmax<=1 excludes evanescent waves)'
  
  write(*,*)'Enter z, the propagation distance'
  read(*,*)z
  write(record_user_inputs_unit,*)z,'  !  z, the propagation distance'
  
  write(*,*)'Do you want to propagate the evansecent field (y/n)?'
  read(*,'(A)')ch
  if ( (ch.eq.'y').OR.(ch.eq.'Y') ) then
    inc_evanescent=.TRUE.
  else
    inc_evanescent=.FALSE.
  end if
  write(record_user_inputs_unit,'(A,A)')ch,'  !  Do you want to propagate the evansecent field (y/n)?'

  w=2.0*pi*frequency
  k=w/c0
  
  pmin=-pmax
  dp=(pmax-pmin)/(np-1)
  
! sort out the numerical integration element of area for 1D and 2D scenarios for the correlation to Wigner function integration
  if (nx.eq.1) then
    dA=2.0*dy
  else if (ny.eq.1) then
    dA=2.0*dx
  else
    dA=2.0*dx*2.0*dy
  end if
      
! sort out the numerical integration element of area for 1D and 2D scenarios for the Wigner function to correlation function integration
  if (nx.eq.1) then
    kdpA=k*dp/(2.0*pi)
  else if (ny.eq.1) then
    kdpA=k*dp/(2.0*pi)
  else
    kdpA=k*k*dp*dp/((2.0*pi)*(2.0*pi))
  end if

  write(*,*)''
  write(*,*)'Wigner function calculation summary'
  write(*,*)''
  write(*,*)'nx=',nx
  write(*,*)'ny=',ny
  write(*,*)'dx=',dx
  write(*,*)'dy=',dy
  write(*,*)'dA=',dA
  write(*,*)'pmin=',pmin
  write(*,*)'pmax=',pmax
  write(*,*)'np=  ',np
  write(*,*)'dp=  ',dp
  write(*,*)'f=   ',frequency,' Hz'
  write(*,*)'k=   ',k,' m^-1'
  write(*,*)'kdpA=',kdpA
  if (inc_evanescent) then
    write(*,*)'Include evanescent field propagation algorithm'
  else
     write(*,*)'No evanescent field propagation algorithm'
  end if
  
! EVALUATE THE WIGNER DISTRIBUTION FUNCTION

! pre-calculate the Fourier integral complex exponentials
  
! may have the range of s wrong here... need to check
  ALLOCATE(exp_function(-ny:ny,-nx:nx,1:np,1:np))

  do ipx=1,npx
  
    if (npx.eq.1) then
      px=0.0
    else
      px=pmin+(ipx-1)*dp
    end if
  
    do ipy=1,npy
  
      if (npy.eq.1) then
    	py=0.0
      else
    	py=pmin+(ipy-1)*dp
      end if
  
      do array_sx=-nx,nx
      	   
    	sx=(real(array_sx))*dx
    	  
  	do array_sy=-ny,ny
  
  	  element_y1=iy1-isy+1
  	  element_y2=isy

  	  sy=(real(array_sy))*dy
    	      
    	  exp_function(array_sy,array_sx,ipy,ipx)=exp(-j*k*(px*sx+py*sy))*dA
    	    
  	end do ! next sy
  
      end do ! next sx
       
    end do ! next py
      
  end do ! next px


  write(*,*)'EVALUATE THE WIGNER DISTRIBUTION FUNCTION'

  Cov2(:,:,:,:)=(0.0,0.0)
      
  write(6,8100,advance='no')'ix1= ',1,' of ',2*nx-1,' iy1= ',1,' of ',2*ny-1
8100    format(A5,I10,A4,I10,A6,I10,A4,I10)    
  flush(6)
  
  do ix1=1,2*nx-1
  
    do iy1=1,2*ny-1
    
      write(6,'(A)',advance='no')char(13)
      write(6,8100,advance='no')'ix1= ',ix1,' of ',2*nx-1,' iy1= ',iy1,' of ',2*ny-1
      flush(6)    
  
      do ipx=1,npx
      
	if (npx.eq.1) then
	  px=0.0
	else
	  px=pmin+(ipx-1)*dp
	end if
  
        do ipy=1,npy
      
	  if (npy.eq.1) then
	    py=0.0
	  else
	    py=pmin+(ipy-1)*dp
	  end if
    
! Integrate over sx
          sum=(0.0,0.0)
    
          do isx=1,ix1
      
            element_x1=ix1-isx+1
            element_x2=isx
	  
! only continue if the element row,col exists : check no 1
	    if (     (element_x1.ge.1).AND.(element_x1.LE.nx) &
	        .AND.(element_x2.ge.1).AND.(element_x2.LE.nx) ) then 
      
! sx2=sx/2 is the speparation in m of the data available in the covariance matrix    	      
              array_sx=element_x2-element_x1
	      
!   Integrate over sy
              do isy=1,iy1
      
                element_y1=iy1-isy+1
                element_y2=isy
	
! only continue if the element row,col exists : check no 2
	        if (     (element_y1.ge.1).AND.(element_y1.LE.ny) &    
	            .AND.(element_y2.ge.1).AND.(element_y2.LE.ny)   ) then 
      
! sy2=sy/2 is the speparation in m of the data available in the covariance matrix  
                  array_sy=element_y2-element_y1
		  
	          sum=sum+exp_function(array_sy,array_sx,ipy,ipx)*Cov(element_x1,element_x2,element_y1,element_y2)
		
	        end if ! check 2 on y element existance
  
              end do ! next sy
	    
	    end if  ! check 1 on x element existance
  
          end do ! next sx
      
          WF(ipy,ipx,iy1,ix1)=sum
      
        end do ! next py
      
      end do ! next px
    
    end do ! next iy
    
  end do ! next ix
  
  DEALLOCATE(exp_function)

! WRITE THE WIGNER FUNCTION TO FILE  

  write(*,*)'WRITE THE WIGNER FUNCTION TO FILE'
  
  write(*,*)'Enter the filename for the output Wigner distribution function data'
  read(*,'(A)')filename
  
  open(unit=local_file_unit,file=trim(filename))
  write(record_user_inputs_unit,'(A)')trim(filename)

  do ix1=1,2*nx-1
  
    x1=sx_values(ix1)
  
    do iy1=1,2*ny-1
    
      y1=sy_values(iy1)
       
      do ipx=1,npx
    
	if (npx.eq.1) then
	  px=0.0
	else
	  px=pmin+(ipx-1)*dp
	end if
      
        do ipy=1,npy
    
	  if (npy.eq.1) then
	    py=0.0
	  else
	    py=pmin+(ipy-1)*dp
	  end if

! If the Wigner function value is too small to fit in the standard Fortran format then set it to zero.
	  if (abs(WF(ipy,ipx,iy1,ix1)).LT.1E-30) WF(ipy,ipx,iy1,ix1)=0.0
      
!          write(local_file_unit,8000)x1,y1,px,py,real(WF(ipy,ipx,iy1,ix1)),Imag(WF(ipy,ipx,iy1,ix1)),Abs(WF(ipy,ipx,iy1,ix1))
          write(local_file_unit,8000)x1,y1,px,py,real(WF(ipy,ipx,iy1,ix1))
8000  format(5E12.4)
      
        end do
      
      end do
    
    end do
    
  end do
  
  close(unit=local_file_unit)

! PROPAGATE THE WIGNER FUNCTION A DISTANCE Z.

  write(*,*)'PROPAGATE THE WIGNER FUNCTION A DISTANCE Z'
  
  WF2(:,:,:,:)=(0.0,0.0)
  Nsamples_WF2(1:np,1:np,1:2*nx-1,1:2*ny-1)=0
  
  do ix1=1,2*nx-1
  
    x1=sx_values(ix1)
  
    do iy1=1,2*ny-1
    
      y1=sy_values(iy1)
       
      do ipx=1,npx
    
	if (npx.eq.1) then
	  px=0.0
	else
	  px=pmin+(ipx-1)*dp
	end if
      
        do ipy=1,npy
    
	  if (npy.eq.1) then
	    py=0.0
	  else
	    py=pmin+(ipy-1)*dp
	  end if
 
          mag_p=sqrt(px*px+py*py)
      
          if (abs(mag_p).LT.1.0) then

! Propagation rule for propagating waves
      
            x2=x1-z*px/sqrt(1.0-abs(mag_p)**2)
            ix2=NINT((x2-xmin)*2.0/dx)
      
            y2=y1-z*py/sqrt(1.0-abs(mag_p)**2)
            iy2=NINT((y2-ymin)*2.0/dy)
      
            if ((ix2.ge.1).AND.(ix2.LE.2*nx-1).AND.(iy2.ge.1).AND.(iy2.LE.2*ny-1)) then
!OLD              WF2(ipy,ipx,iy1,ix1)=WF(ipy,ipx,iy2,ix2)

              WF2(ipy,ipx,iy1,ix1)=WF2(ipy,ipx,iy1,ix1)+WF(ipy,ipx,iy2,ix2)
	      
	      Nsamples_WF2(ipy,ipx,iy1,ix1)=Nsamples_WF2(ipy,ipx,iy1,ix1)+1
	      
            end if
	
          else

! Propagation rule for evanescent waves
! Note that this is still the original algorithm applied here...
            if (inc_evanescent) then
	
              WF2(ipy,ipx,iy1,ix1)=WF(ipy,ipx,iy1,ix1)*exp(-2.0*k*z*sqrt(mag_p*mag_p-1.0))
	  
	    else
	
	      WF2(ipy,ipx,iy1,ix1)=0.0
	
	    end if
      
          end if  ! mag_p.LT.1
      
        end do
    
      end do
      
    end do
    
  end do

! finish the revised algorithm...
  do ix1=1,2*nx-1
  
    do iy1=1,2*ny-1
    
      do ipx=1,npx
          
        do ipy=1,npy
	
	  if(Nsamples_WF2(ipy,ipx,iy1,ix1).NE.0) then
	  
	    WF2(ipy,ipx,iy1,ix1)=WF2(ipy,ipx,iy1,ix1)/Nsamples_WF2(ipy,ipx,iy1,ix1)
	    
	  end if
	
        end do
    
      end do
      
    end do
    
  end do
 
! WRITE THE PROPAGATED WIGNER FUNCTION TO FILE  

  write(*,*)'WRITE THE PROPAGATED WIGNER FUNCTION TO FILE '

  write(*,*)'Enter the filename for the propagated Wigner distribution function data'
  read(*,'(A)')filename
  
  open(unit=local_file_unit,file=trim(filename))
  write(record_user_inputs_unit,'(A)')trim(filename)


  do ix1=1,2*nx-1
  
    x1=sx_values(ix1)
  
    do iy1=1,2*ny-1
    
      y1=sy_values(iy1)
       
      do ipx=1,npx
    
	if (npx.eq.1) then
	  px=0.0
	else
	  px=pmin+(ipx-1)*dp
	end if
      
        do ipy=1,npy
    
	  if (npy.eq.1) then
	    py=0.0
	  else
	    py=pmin+(ipy-1)*dp
	  end if

! If the Wigner function value is too small to fit in the standard Fortran format then set it to zero.
	  if (abs(WF2(ipy,ipx,iy1,ix1)).LT.1D-30) WF2(ipy,ipx,iy1,ix1)=0.0
      
!          write(local_file_unit,8000)x1,y1,px,py,real(WF2(ipy,ipx,iy1,ix1)),Imag(WF2(ipy,ipx,iy1,ix1)),Abs(WF2(ipy,ipx,iy1,ix1))
          write(local_file_unit,8000)x1,y1,px,py,real(WF2(ipy,ipx,iy1,ix1))
      
        end do
      
      end do
    
    end do
    
  end do
  
  close(unit=local_file_unit)
  
! EVALUATE THE COVARIANCE FUNCTION FROM THE WIGNER DISTRIBUTION FUNCTION

  write(*,*)'EVALUATE THE COVARIANCE FUNCTION FROM THE WIGNER DISTRIBUTION FUNCTION'

! pre-calculate the Fourier integral complex exponentials
  ALLOCATE(exp_function(1:np,1:np,-ny:ny,-nx:nx))
    
  do array_sx=-nx,nx
       
    sx=(real(array_sx))*dx
      
    do array_sy=-ny,ny
  
      sy=(real(array_sy))*dy
  
      do ipx=1,npx
     
    	if (npx.eq.1) then
    	  px=0.0
    	else
    	  px=pmin+(ipx-1)*dp
    	end if
  
  	do ipy=1,npy
  
    	  if (npy.eq.1) then
    	    py=0.0
    	  else
    	    py=pmin+(ipy-1)*dp
    	  end if
      
    	  exp_function(ipy,ipx,array_sy,array_sx)=exp(j*k*(px*sx+py*sy))*kdpA
    	  
  	end do ! next py
  
      end do ! next px
    	
    end do ! next sy
  
  end do ! next sx

 
  Cov2(:,:,:,:)=(0.0,0.0)
      
  write(6,8100,advance='no')'ix1= ',1,' of ',2*nx-1,' iy1= ',1,' of ',2*ny-1
  flush(6)

  do ix1=1,2*nx-1
  
    do iy1=1,2*ny-1
    
      write(6,'(A)',advance='no')char(13)
      write(6,8100,advance='no')'ix1= ',ix1,' of ',2*nx-1,' iy1= ',iy1,' of ',2*ny-1
      flush(6)    
    
      do isx=1,ix1   
      
        element_x1=ix1-isx+1
        element_x2=isx
	  
! only continue if the element row,col exists : check no 1
	if (     (element_x1.ge.1).AND.(element_x1.LE.nx) &
	    .AND.(element_x2.ge.1).AND.(element_x2.LE.nx) ) then 
      
! s2=s/2 is the speparation in m of the data available in the covariance matrix  
          array_sx  = element_x2-element_x1
    
          do isy=1,iy1  
      
            element_y1=iy1-isy+1
            element_y2=isy
	
! only continue if the element row,col exists : check no 2
	    if (     (element_y1.ge.1).AND.(element_y1.LE.ny) &    
	        .AND.(element_y2.ge.1).AND.(element_y2.LE.ny)   ) then 
      
! s2=s/2 is the speparation in m of the data available in the covariance matrix    
              array_sy  = element_y2-element_y1
    
! Integrate over px and py
              sum=(0.0,0.0)
  
              do ipx=1,npx
           
                do ipy=1,npy
    	    
	          sum=sum+exp_function(ipy,ipx,array_sy,array_sx)*WF2(ipy,ipx,iy1,ix1)
		  
                end do ! next py
  
              end do ! next px
	    
              Cov2(element_x1,element_x2,element_y1,element_y2)=sum
	  
	    end if ! ! check 2 on y element existance
	    
          end do ! next sy
	  
	end if  ! check 1 on x element existance
        
      end do ! next sx
   
    end do ! next iy
   
  end do ! next ix

  write(*,*)'WRITE THE PROPAGATED COVARIANCE DATA TO FILE'

  write(*,*)'Enter the filename for the output Covariance data'
  read(*,'(A)')filename
  
  open(unit=local_file_unit,file=trim(filename))
  open(unit=local_file_unit2,file=trim(filename)//'_diag')
  write(record_user_inputs_unit,'(A)')trim(filename)

  do ix1=1,nx
    do iy1=1,ny
      do ix2=1,nx
        do iy2=1,ny
    
! If the Covariance function value is too small to fit in the standard Fortran format then set it to zero.
	  if (abs(Cov2(ix1,ix2,iy1,iy2)).LT.1D-30) Cov2(ix1,ix2,iy1,iy2)=0.0
	  
          write(local_file_unit,8010)x_values(ix1),y_values(iy1),x_values(ix2),y_values(iy2),	&
	      real(Cov2(ix1,ix2,iy1,iy2)),Imag(Cov2(ix1,ix2,iy1,iy2)),Abs(Cov2(ix1,ix2,iy1,iy2))
8010  format(7E12.4)

          if ( (ix1.eq.ix2).AND.(iy1.eq.iy2) ) then
! write diagonal elements
            write(local_file_unit2,8020)x_values(ix1),y_values(iy1),	&
	      real(Cov2(ix1,ix2,iy1,iy2)),Imag(Cov2(ix1,ix2,iy1,iy2)),Abs(Cov2(ix1,ix2,iy1,iy2))
	  
8020 format(5E12.4)
	  
	  end if
          
        end do    
      end do   
    end do    
  end do
  
  close(unit=local_file_unit)
  close(unit=local_file_unit2)
 
  DEALLOCATE(Cov)
  DEALLOCATE(WF)
  
  DEALLOCATE(Cov2)
  DEALLOCATE(WF2)
  
  DEALLOCATE( x_values )
  DEALLOCATE( y_values )
  
  DEALLOCATE( sx_values )
  DEALLOCATE( sy_values )
  
  DEALLOCATE(exp_function)


  RETURN
  
END SUBROUTINE wigner_function
