! SUBROUTINE spatial_correlation
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
!    spatial_correlation
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
!     21/1/2014 CJS don't write correlation - no used anywhere really and files are already large...
!
SUBROUTINE spatial_correlation

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character*256	:: ipfilename
  character*256	:: opfilename
  
  integer	:: n_ensemble
  integer	:: n_ensemble_in
  integer	:: nx_points
  integer	:: ny_points
  
  integer	:: n_lines
  
  integer	:: x_col,y_col,re_col,im_col,max_col
  integer	:: step_x,step_y,nxny
  
  integer	:: data_set,xpoint,ypoint,px1,py1,px2,py2
  
  complex*16,allocatable	:: data(:,:,:)
  complex*16,allocatable	:: c(:,:,:,:)
  complex*16,allocatable	:: mean(:,:)
  complex*16			:: cnorm
  
  real*8,allocatable	:: data_line(:)
  
  real*8,allocatable	:: x_value(:)
  real*8,allocatable	:: y_value(:)
  
  real*8	:: x,y,re,im
  real*8	:: x1,y1
  
  logical	:: use_sample_mean
  integer	:: n_samples,i
  
  character	:: ch

! START
  
! READ INPUT FILE

! Need to work out whether we have 1D or 2D data

! Need to work out the number of data sets in the ensemble

  write(*,*)'Enter the filename for the complex ensemble data set'
  read(*,'(A)')ipfilename  
  open(unit=local_file_unit,file=ipfilename)
  write(record_user_inputs_unit,'(A)')trim(ipfilename)

! read the structure of the data file from the user...
  
  write(*,*)'Enter the number of column for x data'
  read(*,*)x_col
  write(record_user_inputs_unit,*)x_col,' Column for x data'
  
  write(*,*)'Enter the number of column for y data'
  read(*,*)y_col
  write(record_user_inputs_unit,*)y_col,' Column for y data'
  
  write(*,*)'Enter the number of column for real part of complex field data'
  read(*,*)re_col
  write(record_user_inputs_unit,*)re_col,' Column for real part of complex field data'
  
  write(*,*)'Enter the number of column for imaginary part of complex field data'
  read(*,*)Im_col
  write(record_user_inputs_unit,*)im_col,' Column for imaginary part of complex field data'
  
  max_col=max(x_col,y_col,re_col,im_col)
  
  ALLOCATE( data_line(1:max_col) )

! Initial read of the data file to work out the size of the scan area in x and y and the number of data sets in the ensemble  
  
  n_lines=0
  nx_points=0
  ny_points=0
  n_ensemble=1
  step_x=0
  step_y=0
  
10 CONTINUE

    read(local_file_unit,*,end=100)(data_line(i),i=1,max_col)
    
    n_lines=n_lines+1
    
    if (n_lines.EQ.1) then
    
! read origin point for spatial scan data
      x1=data_line(x_col)
      y1=data_line(y_col)
      
    else

! work out how the x and y values change in the file
      if ( (data_line(x_col).NE.x1).AND.(step_x.eq.0) ) step_x=n_lines-1
      if ( (data_line(y_col).NE.y1).AND.(step_y.eq.0) ) step_y=n_lines-1 
      
! work out how the number of datasets in the ensemble
      if ( (data_line(x_col).EQ.x1).AND.data_line(y_col).EQ.y1 ) n_ensemble=n_ensemble+1
      
    end if
    
    GOTO 10
    
100 CONTINUE
  
! calculate the number of spatial samples in the x and y directions  

! step_x and step_y should be at least 1...
  if (step_x.eq.0) step_x=n_lines/n_ensemble
  if (step_y.eq.0) step_y=n_lines/n_ensemble

  write(*,*)'Step_x=',step_x
  write(*,*)'Step_y=',step_y
  write(*,*)'n_ensemble=',n_ensemble

! number of points in an xy scan
  nxny=n_lines/n_ensemble
  
  if (step_x.EQ.1) then
    nx_points=step_y
    ny_points=nxny/nx_points
  else if (step_y.eq.1) then
    ny_points=step_x
    nx_points=nxny/ny_points
  else
    write(*,*)'Error in Spatial_correlation: neither step in x and step in y are equal to 1'
    STOP
  end if
  
! do some checks on the file structure  

  write(*,*)'Nx_points =',nx_points
  write(*,*)'Ny_points =',ny_points
  write(*,*)'N_ensemble=',n_ensemble
  write(*,*)'N_lines   =',n_lines

  if (nx_points*ny_points*n_ensemble.NE.n_lines) then
    write(*,*)'Error in Spatial_correlation: nx*ny*n_ensemble.NE.nlines'
    STOP
  end if
    
  write(*,*)'Enter the number of data sets to process or 0 to use the full ensemble'
  read(*,*)n_ensemble_in
  if (n_ensemble_in.NE.0) then
    n_ensemble=n_ensemble_in
  end if
  write(record_user_inputs_unit,*)n_ensemble_in,'  Number of data sets to process or 0 to use the full ensemble'

  write(*,*)'n_ensemble to process=',n_ensemble
  
  ALLOCATE(data(1:nx_points,1:ny_points,1:n_ensemble))
  
  ALLOCATE(mean(1:nx_points,1:ny_points))
 
  ALLOCATE(c(1:nx_points,1:ny_points,1:nx_points,1:ny_points))
  
  ALLOCATE( x_value(1:nx_points) )
  ALLOCATE( y_value(1:ny_points) )
  
  write(*,*)'Do you want to use the sample mean in the correlation calculation (y/n)? (no assumes zero mean)'
  read(*,'(A)')ch
  if ( (ch.EQ.'y').OR.(ch.eq.'Y') ) then
    use_sample_mean=.TRUE.
  else
    use_sample_mean=.FALSE.
  end if
  
  write(record_user_inputs_unit,'(A,A)')ch,'       ! use sample mean in correlation calculation (y/n)'

! Data reading loop  

  rewind(local_file_unit)
  
  do data_set=1,n_ensemble

    if (step_x.eq.1) then
! x point changes fastest
    
      do ypoint=1,ny_points
        do xpoint=1,nx_points
          
          read(local_file_unit,*,end=100)(data_line(i),i=1,max_col)
	  re=data_line(re_col)
	  im=data_line(im_col)
	  x_value(xpoint)=data_line(x_col)
	  y_value(ypoint)=data_line(y_col)
          data(xpoint,ypoint,data_set)=cmplx(re,im)
    
        end do ! next sample in ensemble 
   
      end do ! next ypoint
       
    else
  
      do xpoint=1,nx_points
        do ypoint=1,ny_points       
          
          read(local_file_unit,*,end=100)(data_line(i),i=1,max_col)
	  re=data_line(re_col)
	  im=data_line(im_col)
	  x_value(xpoint)=data_line(x_col)
	  y_value(ypoint)=data_line(y_col)
          data(xpoint,ypoint,data_set)=cmplx(re,im)
    
        end do ! next sample in ensemble 
      end do ! next ypoint
 
    end if ! does x or y increment the fastest
 
  end do ! next dataset in the ensemble

! calculate the mean at each scan point
    
  do xpoint=1,nx_points
  
    do ypoint=1,ny_points
   
      mean(xpoint,ypoint)=(0d0,0d0)
   
      do data_set=1,n_ensemble

        mean(xpoint,ypoint)=mean(xpoint,ypoint)+data(xpoint,ypoint,data_set)
     
      end do
   
      mean(xpoint,ypoint)=mean(xpoint,ypoint)/n_ensemble
   
    end do ! next ypoint
   
  end do ! x point

  if(.NOT.use_sample_mean) then
! *****SET MEAN TO ZERO...******  
    mean(:,:)=(0d0,0d0)
  end if
  
  c(1:nx_points,1:ny_points,1:nx_points,1:ny_points)=(0d0,0d0)
  
  do data_set=1,n_ensemble
    
! calculate the contribution to the correlation matrix from this data set
    do px1=1,nx_points
      do py1=1,ny_points
      
        do px2=1,nx_points
          do py2=1,ny_points
      
            c(px1,py1,px2,py2)=c(px1,py1,px2,py2)+(data(px1,py1,data_set)-mean(px1,py1))*	&
	                                     conjg(data(px2,py2,data_set)-mean(px2,py2))
	
          end do ! next py2
        end do ! next px2
	
      end do ! next py1
    end do ! next px1
    
  end do ! next data set
  
  if(use_sample_mean) then
    n_samples=n_ensemble-1
  else
    n_samples=n_ensemble
  end if
  c(:,:,:,:)=c(:,:,:,:)/(n_samples)

  write(*,*)'Enter the filename for the correlation data'
  read(*,'(A)')opfilename  
  open(unit=12,file=opfilename)
  open(unit=14,file=trim(opfilename)//'_diag')
  write(record_user_inputs_unit,'(A)')trim(opfilename)
  
! write final corelation data
  do px1=1,nx_points
    do py1=1,ny_points
    
      do px2=1,nx_points
  	do py2=1,ny_points
    
! If the Covariance function value is too small to fit in the standard Fortran format then set it to zero.
	  if (abs(c(px1,py1,px2,py2)).LT.1D-30) c(px1,py1,px2,py2)=0d0
   
          cnorm=sqrt(c(px1,py1,px1,py1)*c(px2,py2,px2,py2))
!          write(12,8010)x_value(px1),y_value(py1),x_value(px2),y_value(py2),real(c(px1,py1,px2,py2)),imag(c(px1,py1,px2,py2)),	&
!                        abs(c(px1,py1,px2,py2)),real(c(px1,py1,px2,py2)/cnorm),imag(c(px1,py1,px2,py2)/cnorm)
          write(12,8010)x_value(px1),y_value(py1),x_value(px2),y_value(py2),real(c(px1,py1,px2,py2)),imag(c(px1,py1,px2,py2)),	&
                        abs(c(px1,py1,px2,py2))
8010 format(9E12.4)

          if ( (px1.eq.px2).AND.(py1.eq.py2) ) then
! write diagonal elements
	  
            cnorm=sqrt(c(px1,py1,px1,py1)*c(px2,py2,px2,py2))
!            write(14,8020)x_value(px1),y_value(py1),real(c(px1,py1,px2,py2)),imag(c(px1,py1,px2,py2)),	&
!                        abs(c(px1,py1,px2,py2)),real(c(px1,py1,px2,py2)/cnorm),imag(c(px1,py1,px2,py2)/cnorm)
            write(14,8020)x_value(px1),y_value(py1),real(c(px1,py1,px2,py2)),imag(c(px1,py1,px2,py2)),	&
                        abs(c(px1,py1,px2,py2))
8020 format(7E12.4)
	  
	  end if

  	end do ! next py2
      end do ! next px2

    end do ! next py1
  end do ! next px1
  
  
  close(unit=10)
  close(unit=12)
  close(unit=14)
  
  DEALLOCATE( data_line )
  DEALLOCATE(data)
  DEALLOCATE(mean)
  DEALLOCATE(c)
  DEALLOCATE( x_value )
  DEALLOCATE( y_value )

  RETURN
  
END SUBROUTINE spatial_correlation
