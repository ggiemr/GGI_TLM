! SUBROUTINE Vis_nD
! SUBROUTINE write_4D_vtk_data
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
! SUBROUTINE Vis_nD
!
! NAME
!    Vis_nD
!
! DESCRIPTION
!    create vtk format output file for visualisation of multi-dimensional data 
!    
!     
! COMMENTS
!     We may need to think about an option to interpolate data...
!
! HISTORY
!
!     started 14/08/2014 CJS
!     28/08/2014 CJS: add 1D plot and also make plots 3D
!     28/11/2014 CJS Try to keep the x/y aspect ratio sensible...
!     16/01/2015 CJS Allow the plotting of a reduced set of data (useful for very large data sets)
!     6/02/2015 CJS Change the order of writing points in the vtk file. This improves the look of
!               correlation matrix visualisation especially. 
!     21/1/2016 CJS change to single precision...
!     21/1/2016 CJS allow the code to operate on compressed (.gz) files
!                   also fix a bug with column ranges.
!
SUBROUTINE Vis_nD(function_number)

USE post_process
USE file_information

IMPLICIT NONE

integer	:: function_number

! local variables

  character(len=256)	:: filename
  character(len=256)	:: filename2
  
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
  integer		:: n_columns,column,data_column
  integer		:: data_col
  
  integer		:: loop,col_loop
  
  integer		:: sample,n_samples
  
  integer,allocatable	:: column_list(:)
  
  real,allocatable	:: data(:,:)
  real,allocatable	:: max_data(:)
  real,allocatable	:: min_data(:)
  
  real		:: min_value,value_range,zrange
  
  real		:: value(4)
  real		:: new_value(4)
  integer		:: step(0:4)
  integer		:: n(0:4)

  real,allocatable    :: x_values(:)
  real,allocatable    :: y_values(:)

  real,allocatable    :: plot_data(:,:)
 
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
  
  logical :: write_map

! compression stuff
  logical	:: compression_flag
  integer	:: len_filename
  character(len=256)	:: gzfilename
  character*3   :: extn

! START

  DATA axis_label / 'x_small','y_small','x_large','y_large' /  

  write(*,*)
  write(*,*)'Output files:'
  command='ls -ltr '
  CALL system(command)

! get the filename for the data
5 CONTINUE
  write(*,*)'Enter the filename for the input data'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)

! check for .gz extension which indicates a compressed file
  compression_flag=.FALSE.
  len_filename=LEN(trim(filename))
  extn=filename(len_filename-2:len_filename)
  if (extn.EQ.'.gz') then
    compression_flag=.TRUE.
  end if
  
! open and read the file
  
  if (compression_flag) then
! We must give the filename without the .gz extension here
    gzfilename=filename(1:len_filename-3)
    CALL open_output_file_read(local_file_unit,gzfilename,compression_flag)
  else
    OPEN(unit=local_file_unit,file=filename)
  end if
  
  CALL write_file_format_information(local_file_unit,n_lines,number_of_comment_lines,number_of_data_lines)
  
  rewind(unit=local_file_unit)
  
  write(*,*)'Enter the first line of the numeric data to plot'
  read(*,*)first_line
  write(record_user_inputs_unit,*)first_line,' First line of the numeric data to plot'
  
  write(*,*)'Enter the last line of the data file to process or 0 to read the whole file'
  read(*,*)last_line
  write(record_user_inputs_unit,*)last_line,' Last line of the numeric data to plot or 0 to read the whole file'
  
  if ( (last_line.EQ.0).OR.(last_line.GT.(n_lines-number_of_comment_lines)) ) then
    last_line=n_lines-number_of_comment_lines
  end if 
  
  if (last_line.LT.first_line) then
    write(*,*)'Error: last line specified is less than the first line'
    STOP
  end if 
      
  write(*,*)'Enter the total number of columns of data to read (including coordinate data and data value)'
  read(*,*)max_columns
  write(record_user_inputs_unit,*)max_columns,' Number of columns of data to read (including coordinate data and data value)'
    
  n_samples=last_line-first_line+1
    
  ALLOCATE ( data(1:n_samples,1:max_columns) )
  ALLOCATE ( max_data(1:max_columns) )
  ALLOCATE ( min_data(1:max_columns) )
      
  max_data(1:max_columns)=-1D30
  min_data(1:max_columns)= 1D30

! read lines to ignore
  do i=1,number_of_comment_lines
     read(local_file_unit,*,ERR=9000)
  end do
    
  do sample=1,n_samples
    
    read(local_file_unit,*,ERR=9000)(data(sample,i),i=1,max_columns)

    do loop=1,max_columns
	  
      max_data(loop)=max( max_data(loop),data(sample,loop) )
      min_data(loop)=min( min_data(loop),data(sample,loop) )
	  
    end do
     
  end do ! next sample to read

  if (compression_flag) then
    CALL close_output_file(local_file_unit,filename2,compression_flag)
  else
    CLOSE(unit=local_file_unit)
  end if
 
! Write column data ranges to the screen      
  write(*,*)'Number of data samples read:',n_samples
  write(*,*)'Data column ranges:'
  do col_loop=1,max_columns
    write(*,'(A,I3,A,E16.6,A,E16.6)')'Column ',col_loop,' min value=',min_data(col_loop),' max value=',max_data(col_loop)
  end do
  
! organise the data into a sensible format
  write(*,*)"Enter the number of dimensions to plot (between 1 and 4) and add 'plus_map' to plot the tile outlines"
  write(*,*)"e.g. to generate a 4D plot with tile outlines enter:"
  write(*,*)"4 plus_map"
  
  line=''
  read(*,'(A)')line
  read(line,*)nDim
  
  if (index(line,'plus_map').NE.0) then
    write_map=.TRUE.
    write(record_user_inputs_unit,*)nDim,'plus_map : Number of dimensions to plot (between 1 and 4)'
  else
    write_map=.FALSE.
    write(record_user_inputs_unit,*)nDim,' Number of dimensions to plot (between 1 and 4)'
  end if
      
  ALLOCATE( column_list(1:nDim+1) )
  
  data_column=nDim+1
  
  do loop=1,nDim
  
    write(*,'(A,A,A)')'Which Column should the ',axis_label(loop),' axis data come from? '
    read(*,*)column_list(loop)
    write(record_user_inputs_unit,*)column_list(loop),' Column for the ',axis_label(loop),' axis data '
    
  end do 
  
  write(*,'(A,A,A)')'Which Column should the plot value data come from? '
  read(*,*)column_list(data_column)
  write(record_user_inputs_unit,*)column_list(data_column),' Column for the plot value data '
  
  min_value=min_data(column_list(data_column))
  value_range=(max_data(column_list(data_column))-min_data(column_list(data_column)))
   
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
  
  do i=1,nDim
  
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
  
  do i=1,nDim
    write(*,'(A,I1,A,A,A,I10)')'n(',i,')=n(',axis_label(i),')=',n(i)
  end do
      
! get the plotting dimensions
  nx=n(1)*n(3)
  ny=n(2)*n(4)
  
! ny=1 for 1D plots, we must increase this in order to generate a surface to plot
  if (ndim.eq.1) ny=2
  
  write(*,*)'nx=',nx
  write(*,*)'ny=',ny
   
  if (ndim.GT.2) then 
  
    write(*,*)'Do you want to restrict the range of the plot (i.e. number of sub plots) (y/n)?'
    read(*,'(A)')ch
    write(record_user_inputs_unit,'(A,A)')ch,' ! Restrict range of plot (i.e. number of sub plots)'
    
  else
  
! only a 1D or 2D plot so don't restrict the number of plots  
    ch='n'
    
  end if
  
  if ( (ch.EQ.'y').OR.(ch.EQ.'Y') ) then

! Only restrict the large x and large y axis ranges
    nplotmin(1:2)=1
    nplotmax(1:2)=n(1:2)
      
! Large x range
    write(*,*)'Number of plots in the large_x range is',n(3)
    
    write(*,*)'Enter the first plot of the xrange to include'
    read(*,*)nplotmin(3)
    write(record_user_inputs_unit,*)nplotmin(3),' First plot of the xrange to include'
    
    write(*,*)'Enter the last plot of the xrange to include'
    read(*,*)nplotmax(3)
    write(record_user_inputs_unit,*)nplotmax(3),' Last plot of the xrange to include'
      
! Large y range
    write(*,*)'Number of plots in the large_y range is',n(4)
    
    write(*,*)'Enter the first plot of the yrange to include'
    read(*,*)nplotmin(4)
    write(record_user_inputs_unit,*)nplotmin(4),' First plot of the yrange to include'
    
    write(*,*)'Enter the last plot of the yrange to include'
    read(*,*)nplotmax(4)
    write(record_user_inputs_unit,*)nplotmax(4),' Last plot of the yrange to include'

! get the new plotting dimensions
    nxplot=(nplotmax(1)-nplotmin(1)+1)*(nplotmax(3)-nplotmin(3)+1)
    nyplot=(nplotmax(2)-nplotmin(2)+1)*(nplotmax(4)-nplotmin(4)+1)
   
  else
! don't restrict plot range
    nplotmin(1:4)=1
    nplotmax(1:4)=n(1:4)
    nxplot=nx
    nyplot=ny
  end if

! get the new plotting dimensions
 
  ALLOCATE( x_values(1:nxplot) )
  ALLOCATE( y_values(1:nyplot) )
  ALLOCATE( plot_data(1:nxplot,1:nyplot) )

! define the plotting array - this is rather aritrary at the moment: 
! keep the x/y aspect ratio unless it is too large (too small)
  
  if ( ( (nxplot*10.LT.nyplot).OR.(nyplot*10.LT.nxplot) ).AND.(ndim.ne.1) ) then

! large aspect ratio- make the plot square
    dx=1d0/(nxplot-1)
    dy=1d0/(nyplot-1)
    write(*,*)'Large aspect ratio, make plot square'
    
  else
! aspect ratio less than 1:5 or 1D plot so keep it as it is

    dx=1d0
    dy=1d0
    
  end if
  
  lx=dx*(nx-1)
  ly=dy*(ny-1)
  
! set the x and y values for the surface plot
  write(*,*)'nxplot=',nxplot,' dx=',dx,' lx=',lx
  write(*,*)'nyplot=',nyplot,' dy=',dy,' ly=',ly

! loop over the whole x axis and save the x data only where we are plotting results
  ix=0
  i=0
  do a3=1,n(3)  ! large x axis
    do a1=1,n(1)  ! small x axis
    
      ix=ix+1
      
      if ( (a3.GE.nplotmin(3)).AND.(a3.LE.nplotmax(3)) ) then
        i=i+1
        x_values(i)=(ix-1)*dx
      end if
      
    end do
  end do
  
! loop over the whole y axis and save the y data only where we are plotting results
  iy=0
  i=0
  do a4=1,n(4)  ! large y axis
    do a2=1,n(2)  ! small y axis
    
      iy=iy+1
      if ( (a4.GE.nplotmin(4)).AND.(a4.LE.nplotmax(4)) ) then
        i=i+1
        y_values(i)=(iy-1)*dy
      end if
      
    end do
  end do
 
  write(*,*)'Get the plot data'

  iy=0
  
  do a4=nplotmin(4),nplotmax(4)  ! large y axis
    do a2=nplotmin(2),nplotmax(2)  ! small y axis

      iy=iy+1
      
      ix=0
      
      do a3=nplotmin(3),nplotmax(3)  ! large x axis
        do a1=nplotmin(1),nplotmax(1)  ! small x axis
      
          ix=ix+1
	  
	  element(4)=a4
	  element(3)=a3
	  element(2)=a2
	  element(1)=a1

! checks on the array ranges	  
	  do i=1,4
	    if (element(i).GT.n(i)) then
	      write(*,*)'Dimension ',i,' error'
	      write(*,*)'element(i)=',element(i),' n(i)=',n(i)
	      STOP
	    end if
          end do
	  
	  if (ix.GT.nxplot) then
	    write(*,*)'ix error',ix,nxplot
	    STOP
	  end if
	  
	  if (iy.GT.nyplot) then
	    write(*,*)'iy error',iy,nyplot
	    STOP
	  end if
	  
	  sample=1+(a1-1)*step(1)+(a2-1)*step(2)+(a3-1)*step(3)+(a4-1)*step(4)
	  
	  if (sample.GT.n_samples) then
	  
	    write(*,*)'sample error',sample,n_samples
	    write(*,*)a1,n(1),step(1)
	    write(*,*)a2,n(2),step(2)
	    write(*,*)a3,n(3),step(3)
	    write(*,*)a4,n(4),step(4)
	    STOP
	  end if
	  
          plot_data(ix,iy)=Data(sample,column_list(data_column))
	  	  
        end do    
      end do   
    end do    
  end do
  
! again a special case for 1 dimensional data, artificially add another column of y data
  if (ndim.eq.1) then
    plot_data(1:nxplot,2)=plot_data(1:nxplot,1)
  end if

! set the height (zrange) of the surface plot to be related to the x and y extent of the plot  
  zrange=max(dx*(nxplot-1),dy*(nyplot-1))/4d0
  
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

  filename2=trim(filename)//".vtk"
  CALL write_4D_vtk_data(nxplot,nyplot,x_values,y_values,plot_data,min_value,value_range,zrange,filename2)    
  
  DEALLOCATE( x_values )
  DEALLOCATE( y_values )
  DEALLOCATE( plot_data )

! write a file which can be used to show the edges of the small plots
  
  nx=n(3)+1
  ny=n(4)+1

  ALLOCATE( x_values(1:nx) )
  ALLOCATE( y_values(1:ny) )
  ALLOCATE( plot_data(1:nx,1:ny) )
  
  if (ndim.ne.1) then
  
    x_values(1)=-0.5*dx
    do i=2,nx
      x_values(i)=x_values(i-1)+n(1)*dx
    end do

    y_values(1)=-0.5*dy
    do i=2,ny
      y_values(i)=y_values(i-1)+n(2)*dy
    end do
   
  else if (ndim.eq.1) then
! Special case for 1 dimensional data
     
    x_values(1)=-0.5
    do i=2,nx
      x_values(i)=x_values(i-1)+n(1)
    end do

    y_values(1)=-0.5
    y_values(2)=1.5
   
  end if
  
  if (write_map) then
    plot_data(1:nx,1:ny)=0d0
    min_value=0d0
    value_range=1d0
    zrange=1d0
    filename2=trim(filename)//".map.vtk"
    CALL write_4D_vtk_data(nx,ny,x_values,y_values,plot_data,min_value,value_range,zrange,filename2)    
  end if
   
  DEALLOCATE ( data )
  DEALLOCATE ( max_data )
  DEALLOCATE ( min_data )
  DEALLOCATE( column_list )
  DEALLOCATE( x_values )
  DEALLOCATE( y_values )
  DEALLOCATE( plot_data )

  RETURN
     
9000 CALL write_line('Error reading data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(filename)
     STOP 
  
END SUBROUTINE Vis_nD
!
! _______________________________________________________________
!
!    
SUBROUTINE write_4D_vtk_data(nx,ny,x_values,y_values,data,min_value,value_range,zrange,filename)    

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

END SUBROUTINE write_4D_vtk_data
