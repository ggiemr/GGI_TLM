! SUBROUTINE animate_near_field_time_data
! SUBROUTINE write_TD_vtk_data
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
! SUBROUTINE animate_near_field_time_data
!
! NAME
!    animate_near_field_time_data
!
! DESCRIPTION
!    create vtk format output files for animation time domain near field scan data
!    
!     
! COMMENTS
!    
!
! HISTORY
!
!     started 6/11/2024 CJS
!
SUBROUTINE animate_near_field_time_data(function_number)

USE constants
USE post_process
USE file_information

IMPLICIT NONE

integer	:: function_number

! local variables

  character(len=256)	:: filename
  character(len=256)	:: filename2
  character(len=256)	:: frame_filename
  
  character(len=256)	:: command
  character(len=256)	:: line
  
  logical		:: file_exists
	   
  integer   		:: first_line
  integer   		:: last_line
  integer   		:: n_lines
  
  integer :: number_of_comment_lines
  integer :: number_of_data_lines

  integer		:: nDim
  integer		:: Dim(4)
  
  character*7		:: axis_label(4)
  
  integer		:: max_columns
  integer		:: n_columns,column,Re_data_column,Im_data_column
  integer		:: data_col
  
  integer		:: loop,col_loop
  
  integer		:: sample,n_samples
  
  integer,allocatable	:: column_list(:)
  
  real,allocatable	:: data(:,:)
  real,allocatable	:: data_line(:)
  real,allocatable	:: max_data(:)
  real,allocatable	:: min_data(:)
  
  real		:: min_value,value_range
  real		:: min_Re_value,Re_value_range
  real		:: min_Im_value,Im_value_range
  real          :: zrange
  
  real		:: value(4)
  real		:: new_value(4)
  integer		:: step(0:4)
  integer		:: n(0:4)

  real,allocatable    :: x_values(:)
  real,allocatable    :: y_values(:)

  real,allocatable    :: plot_data(:,:)
  complex,allocatable    :: cmplx_plot_data(:,:)
 
  integer		:: i1,i2,i3,i4
  integer		:: a1,a2,a3,a4

  integer               ::  nx,ny
  integer               ::  ix,iy
  
  integer		:: nxplot,nyplot
  integer		:: nplotmin(4),nplotmax(4)

  integer		:: element(4)
 
  integer		:: i
  real		:: dx,dy,lx,ly
  
  integer 		:: num,den
    
  character		:: ch
  character(LEN=256)    :: file_line
  
  logical :: write_map

! compression stuff
  logical	:: compression_flag
  integer	:: len_filename
  character(len=256)	:: gzfilename
  character*3   :: extn
  
! animation stuff

  integer :: frame,n_frames
  
  real*8 :: re,im,phase

! START

  DATA axis_label / 'x_small','y_small','x_large','y_large' /  
  
  write(*,*)
  write(*,*)'Output files:'
  command='ls -ltr '
  CALL system(command)

! get the filename for the data
5 CONTINUE
  write(*,*)'Enter the filename for the near field time domain data'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)

  OPEN(unit=local_file_unit,file=filename)

! the format is known: x y z t Vx Vy Vz  so read 7 columns of data

  CALL write_file_format_information(local_file_unit,n_lines,number_of_comment_lines,number_of_data_lines)
    
  rewind(unit=local_file_unit)
  
  first_line=1    ! 
  
  last_line=0     ! read to the last line
  
  if ( (last_line.EQ.0).OR.(last_line.GT.(n_lines-number_of_comment_lines)) ) then
    last_line=n_lines-number_of_comment_lines
  end if 
  
  if (last_line.LT.first_line) then
    write(*,*)'Error: last line specified is less than the first line'
    STOP
  end if 
  
  write(*,*)'First line to read=',first_line
  write(*,*)'Last line to read=',last_line
      
  max_columns=7
  write(*,*)max_columns,' Number of columns of data to read (including coordinate data and data values)'
    
  n_samples=last_line-first_line+1
  write(*,*)'Number of samples to read=',n_samples
  
  write(*,*)'Allocate memory for the data'
  ALLOCATE ( data_line(1:max_columns) )
  ALLOCATE ( data(1:n_samples,1:max_columns) )
  ALLOCATE ( max_data(1:max_columns) )
  ALLOCATE ( min_data(1:max_columns) )
      
  max_data(1:max_columns)=-1D30
  min_data(1:max_columns)= 1D30

! read the initial line to ignore
  read(local_file_unit,*,ERR=9000)
  
  sample=0

10  CONTINUE

  read(local_file_unit,'(A256)',END=100)file_line
        
  read(file_line,*,ERR=10,END=10)(data_line(i),i=1,max_columns)  ! read the data - read the next line if there is an error

  sample=sample+1
  
  if (sample.LE.n_samples) then 
    do loop=1,max_columns
      data(sample,loop)=data_line(loop)
      max_data(loop)=max( max_data(loop),data(sample,loop) )
      min_data(loop)=min( min_data(loop),data(sample,loop) )
	  
    end do
  end if
  
GOTO 10  ! read the next line

100 CONTINUE

CLOSE(unit=local_file_unit)

write(*,*)'Number of samples read    =',sample
write(*,*)'Number of samples expected=',n_samples
  
 
! Write column data ranges to the screen      
  write(*,*)'Number of data samples read:',n_samples
  write(*,*)'Data column ranges:'
  do col_loop=1,max_columns
    write(*,'(A,I3,A,E16.6,A,E16.6)')'Column ',col_loop,' min value=',min_data(col_loop),' max value=',max_data(col_loop)
  end do
  
  n_frames=NINT(max_data(4))-NINT(min_data(4))+1
  write(*,*)'Number of time frames=',n_frames
  
! organise the data into a sensible format

  nDim=3
      
  ALLOCATE( column_list(1:nDim+1) )
  
  Re_data_column=nDim+1
    
  write(*,'(A,A,A)')'Which Column should the x axis data come from? '
  read(*,*)column_list(1)
  write(record_user_inputs_unit,*)column_list(1),' Column for the x axis data '
    
  write(*,'(A,A,A)')'Which Column should the y axis data come from? '
  read(*,*)column_list(2)
  write(record_user_inputs_unit,*)column_list(2),' Column for the y axis data '
    
  write(*,'(A,A,A)')'Which Column should the time axis data come from? '
  read(*,*)column_list(3)
  write(record_user_inputs_unit,*)column_list(3),' Column for the time axis data '
    
  write(*,'(A,A,A)')'Which Column should the Real plot value data come from? '
  read(*,*)column_list(Re_data_column)
  write(record_user_inputs_unit,*)column_list(Re_data_column),' Column for the Real plot value data '
    
  min_Re_value=min_data(column_list(Re_data_column))
  Re_value_range=(max_data(column_list(Re_data_column))-min_data(column_list(Re_data_column)))
   
! work out the size of each dimension in the file

  step(:)=0

! get initial plot axis values
  value(:)=0d0
  
  do loop=1,nDim
    value(loop)=data(1,column_list(loop)) 
  end do

! keep reading through the file looking for changes in each of the dimension variables
  do sample=2,n_samples
  
    do loop=1,nDim
      new_value(loop)=data(sample,column_list(loop)) 
    end do
  
    do i=1,nDim
      if ( (new_value(i).NE.value(i)).AND.(step(i).eq.0) ) step(i)=sample-1
    end do

  end do
  
  step(0)=n_samples
  
! work out the number of elements in each dimensions
  
  do i=0,nDim
    write(*,'(A,I2,A,I10)')'Step ',i,'=',step(i)
  end do
  
  n(:)=1 ! default dimension is 1
  
  do i=1,2
  
    den=step(i)
! find the next largest step value
    num=999999999
    do loop=0,nDim
      if ( (step(loop).GT.den).AND.(step(loop).LT.num) ) then
        num=step(loop)
      end if
    end do
    n(i)=num/den
  end do
  
  do i=1,2
    write(*,'(A,I1,A,A,A,I10)')'n(',i,')=n(',axis_label(i),')=',n(i)
  end do
      
! get the plotting dimensions
  nx=n(1)
  ny=n(2)
    
  write(*,*)'nx=',nx
  write(*,*)'ny=',ny
   
! don't restrict plot range

  nplotmin(1:2)=1
  nplotmax(1:2)=n(1:2)
  nxplot=nx
  nyplot=ny

! get the plotting dimensions
 
  ALLOCATE( x_values(1:nxplot) )
  ALLOCATE( y_values(1:nyplot) )
  ALLOCATE( plot_data(1:nxplot,1:nyplot) )

  sample=1
  do ix=1,nx
    x_values(ix)=data(sample,column_list(1))
    sample=sample+1
    write(*,*)'ix=',ix,' x=',x_values(ix)
  end do

  sample=1
  do iy=1,ny
    y_values(iy)=data(sample,column_list(2))
    sample=sample+nx
    write(*,*)'iy=',iy,' y=',y_values(iy)
  end do

  dx=(x_values(nx)-x_values(1))/real(nx-1)
  dy=(y_values(ny)-y_values(1))/real(ny-1)
  
  lx=dx*(nx-1)
  ly=dy*(ny-1)
  
! set the x and y values for the surface plot
  write(*,*)'nxplot=',nxplot,' dx=',dx,' lx=',lx
  write(*,*)'nyplot=',nyplot,' dy=',dy,' ly=',ly
  
! set the height (zrange) of the surface plot to be related to the x and y extent of the plot  
  zrange=max(dx*(nxplot-1),dy*(nyplot-1))/4d0
  
  min_value=min_Re_value
  value_range=Re_value_range
  
  write(*,*)'minimum Real data value=',min_Re_value
  write(*,*)'Real value_range       =',Re_value_range
  write(*,*)'minimum Imaginary data value=',min_Im_value
  write(*,*)'Imaginary value_range       =',Im_value_range
  write(*,*)'minimum data value=',min_value
  write(*,*)'value_range       =',value_range
  write(*,*)'z_range           =',zrange
 
  write(*,*)'Enter the filename visualisation data (without vtk extension) or change_range'
  read(*,'(A)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  if (filename(1:12).EQ.'change_range') then
  
! Alter the range of the plot

    write(*,*)'Enter the minimum value for the data plot'
    read(*,*)min_value
    write(record_user_inputs_unit,*)min_value
    
    write(*,*)'Enter the value range for the data plot'
    read(*,*)value_range
    write(record_user_inputs_unit,*)value_range
    
    write(*,*)'Enter the z range for the data plot'
    read(*,*)zrange
    write(record_user_inputs_unit,*)zrange
    
    write(*,*)'Enter the filename visualisation data (without vtk extension)'
    read(*,'(A)')filename
    write(record_user_inputs_unit,'(A)')trim(filename)
    
  end if

! START ANIMATION LOOP
  
  sample=0
  
  do frame=1,n_frames
      
    write(*,*)'Writing frame data',frame
    
! add frame number to base filename

    CALL add_integer_to_filename(filename,frame,filename2)

    frame_filename=trim(filename2)//'.vtk'
    
    write(*,*)'Writing frame to file:',trim(frame_filename)
  
! get the values to plot 

    do iy=1,ny
      do ix=1,nx
        sample=sample+1
        plot_data(ix,iy)=data(sample,column_list(Re_data_column))
      end do
    end do
 
! Write the data to vtk file

    CALL write_TD_vtk_data(nxplot,nyplot,x_values,y_values,plot_data,min_value,value_range,zrange,frame_filename)    
  
  end do ! next frame 
  
  DEALLOCATE( x_values )
  DEALLOCATE( y_values )
  DEALLOCATE( plot_data )

   
  DEALLOCATE ( data )
  DEALLOCATE ( data_line )
  DEALLOCATE ( max_data )
  DEALLOCATE ( min_data )
  DEALLOCATE( column_list )

  RETURN
     
9000 CALL write_line('Error reading data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(filename)
     STOP 
  
END SUBROUTINE animate_near_field_time_data
!
! _______________________________________________________________
!
!    
SUBROUTINE write_TD_vtk_data(nx,ny,x_values,y_values,data,min_value,value_range,zrange,filename)    

  IMPLICIT NONE

  integer	::  nx,ny
  
  real	:: x_values(1:nx)
  real	:: y_values(1:ny)
  
  real	:: data(1:nx,1:ny)
  
  real	:: min_value,value_range,zrange
  
  character*256 :: filename

! local variables

  integer	:: n_points,point
  integer	:: n_quads,quad
  
  integer	:: ix,iy
  
  integer	:: count

! START
    
  open(UNIT=20,FILE=trim(filename))
  
  n_points=nx*ny
  n_quads=(nx-1)*(ny-1)
  
  write(*,*)'Writing vtk data'
  write(*,*)'nx      =',nx
  write(*,*)'ny      =',ny
  write(*,*)'n_points=',n_points
  write(*,*)'n_quads =',n_quads
    
! write header information    
! write vtk header      
  write(20,'(A)')'# vtk DataFile Version 2.0'
  write(20,'(A)')trim(filename)
  write(20,'(A)')'ASCII'
  write(20,'(A)')'DATASET POLYDATA'
  write(20,'(A,I10,A)')'POINTS',n_points,' float'
    
! write point data 
    
  count=0
    
  do ix=1,nx
    do iy=1,ny
      
      write(20,8000)x_values(ix),y_values(iy),(data(ix,iy)-min_value)*zrange/value_range
8000  format(3E12.4)
      count=count+1
      
    end do ! next iy
  end do ! next ix
  
  write(*,*)'Number of points written=',count,' n_points=',n_points
  
! write quad data

  write(*,*)'Write quad data'
  write(20,'(A,2I10)')'POLYGONS',n_quads,n_quads*5

  point=0
  count=0
  
  do ix=1,nx-1
  
    do iy=1,ny-1
      
!      write(20,8010)4,point,point+1,point+ny+1,point+ny
! this ordering of the points may lead to a better interpolation of some datasets...
      write(20,8010)4,point+1,point,point+ny,point+ny+1
8010  format(I3,4I8)
      count=count+1

      point=point+1

    end do ! next iy
    
    point=point+1
     
  end do ! next ix
  
  write(*,*)'Number of quads written=',count,' n_quads=',n_quads
  
! STAGE 6. write data associated with points

! write point based data

  write(*,*)'Write point data'

  write(20,'(A,I10)')'POINT_DATA ',n_points
  write(20,'(A)')'SCALARS Field_on_cells float 1'
  write(20,'(A)')'LOOKUP_TABLE field_on_cells_table'
  
  count=0
  do ix=1,nx
  
    do iy=1,ny
      
      write(20,8020)data(ix  ,iy  )
      count=count+1
      
8020  format(E12.4)

    end do ! next iy
     
  end do ! next ix
  
  write(*,*)'Number of data points written=',count,' n_points=',n_points

  close(UNIT=20)

END SUBROUTINE write_TD_vtk_data
