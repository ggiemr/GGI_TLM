! SUBROUTINE animate_2D_complex_data
! SUBROUTINE write_2D_vtk_data
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
! SUBROUTINE animate_2D_complex_data
!
! NAME
!    animate_2D_complex_data
!
! DESCRIPTION
!    create vtk format output files for animation of 2D complex data 
!    the data may be on a non-rectangular grid
!    
!     
! COMMENTS
!    based on animate_nD_complex_data
!
! HISTORY
!
!     started 22/1/2025 CJS
!
SUBROUTINE animate_2D_complex_data(function_number)

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
  integer		:: Dim(2)
  
  character*7		:: axis_label(2)
  
  integer		:: max_columns
  integer		:: n_columns,column
  integer		:: i_data_column,j_data_column
  integer		:: x_data_column,y_data_column
  integer		:: Re_data_column,Im_data_column
  integer		:: data_col
  
  integer		:: loop,col_loop
  
  integer		:: sample,n_samples
  
  integer,allocatable	:: column_list(:)
  
  real,allocatable	:: data(:,:)
  real,allocatable	:: max_data(:)
  real,allocatable	:: min_data(:)
  
  real		:: min_value,value_range
  real		:: min_Re_value,Re_value_range
  real		:: min_Im_value,Im_value_range
  real          :: zrange
  
  real		:: value(2)
  real		:: new_value(2)
  integer		:: step(0:2)
  integer		:: n(0:2)

  real,allocatable    :: i_values(:)
  real,allocatable    :: j_values(:)
  real,allocatable    :: xy_values(:,:,:)
  
  real :: max_xy,min_xy

  real,allocatable    :: plot_data(:,:)
  complex,allocatable    :: cmplx_plot_data(:,:)
 
  integer		:: i1,i2
  integer		:: a1,a2

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
  
! animation stuff

  integer :: frame,n_frames
  
  real*8 :: re,im,phase

! START

  DATA axis_label / 'x','y' /  
  
  n_frames=20

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
      
  write(*,*)'Enter the total number of columns of data to read (including coordinate data and data values)'
  read(*,*)max_columns
  write(record_user_inputs_unit,*)max_columns,' Number of columns of data to read (including coordinate data and data values)'
    
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
  
  nDim=2
  
  write_map=.FALSE.
      
  ALLOCATE( column_list(1:6) )
  
  i_data_column=1
  j_data_column=2
  x_data_column=3
  y_data_column=4
  Re_data_column=5
  Im_data_column=6
    
  write(*,'(A)')'Which Column should the i axis data come from? '
  read(*,*)column_list(i_data_column)
  write(record_user_inputs_unit,*)column_list(i_data_column),' Column for the i axis data '
    
  write(*,'(A)')'Which Column should the j axis data come from? '
  read(*,*)column_list(j_data_column)
  write(record_user_inputs_unit,*)column_list(j_data_column),' Column for the j axis data '
    
  write(*,'(A)')'Which Column should the x axis data come from? '
  read(*,*)column_list(x_data_column)
  write(record_user_inputs_unit,*)column_list(x_data_column),' Column for the x axis data '
    
  write(*,'(A)')'Which Column should the y axis data come from? '
  read(*,*)column_list(y_data_column)
  write(record_user_inputs_unit,*)column_list(y_data_column),' Column for the y axis data '
  
  write(*,'(A)')'Which Column should the Real plot value data come from? '
  read(*,*)column_list(Re_data_column)
  write(record_user_inputs_unit,*)column_list(Re_data_column),' Column for the Real plot value data '
  
  write(*,'(A)')'Which Column should the Imaginary plot value data come from? '
  read(*,*)column_list(Im_data_column)
  write(record_user_inputs_unit,*)column_list(Im_data_column),' Column for the Imaginary plot value data '
  
  min_Re_value=min_data(column_list(Re_data_column))
  Re_value_range=(max_data(column_list(Re_data_column))-min_data(column_list(Re_data_column)))
  
  min_Im_value=min_data(column_list(Im_data_column))
  Im_value_range=(max_data(column_list(Im_data_column))-min_data(column_list(Im_data_column)))
   
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
  nx=n(1)
  ny=n(2)
    
  write(*,*)'nx=',nx
  write(*,*)'ny=',ny
   
  nplotmin(1:2)=1
  nplotmax(1:2)=n(1:2)
  nxplot=nx
  nyplot=ny

! get the new plotting dimensions
 
  ALLOCATE( i_values(1:nxplot) )
  ALLOCATE( j_values(1:nyplot) )
  ALLOCATE( xy_values(1:2,1:nxplot,1:nyplot) )
  ALLOCATE( cmplx_plot_data(1:nxplot,1:nyplot) )
  ALLOCATE( plot_data(1:nxplot,1:nyplot) )

! define the plotting array - this is rather aritrary at the moment: 
! keep the x/y aspect ratio unless it is too large (too small)
  
  dx=1d0
  dy=1d0
  
  lx=dx*(nx-1)
  ly=dy*(ny-1)
  
! set the x and y values for the surface plot
  write(*,*)'nxplot=',nxplot,' dx=',dx,' lx=',lx
  write(*,*)'nyplot=',nyplot,' dy=',dy,' ly=',ly

! loop over the whole x axis and save the x data only where we are plotting results
  ix=0
  do a1=1,n(1)  ! x axis
    
    ix=ix+1  
    i_values(ix)=(ix-1)*dx
      
  end do
  
! loop over the whole y axis and save the y data only where we are plotting results
  iy=0
  do a2=1,n(2)  ! small y axis
    
    iy=iy+1
    j_values(i)=(iy-1)*dy
      
  end do
 
  write(*,*)'Get the complex plot data'
  
  iy=0
  
  do a2=nplotmin(2),nplotmax(2)  ! y axis

    iy=iy+1
    
    ix=0
    
    do a1=nplotmin(1),nplotmax(1)  ! small x axis
    
      ix=ix+1
        
      element(2)=a2
      element(1)=a1

! checks on the array ranges	  
      do i=1,2
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

      sample=1+(a1-1)*step(1)+(a2-1)*step(2)

      if (sample.GT.n_samples) then
        write(*,*)'sample error',sample,n_samples
        write(*,*)a1,n(1),step(1)
        write(*,*)a2,n(2),step(2)
        STOP
      end if
      
      re=Data(sample,column_list(Re_data_column))
      im=Data(sample,column_list(Im_data_column))
          
      cmplx_plot_data(ix,iy)=dcmplx(re,im)
      xy_values(1,ix,iy)=Data(sample,column_list(x_data_column))
      xy_values(2,ix,iy)=Data(sample,column_list(y_data_column))
          
    end do    
  end do
  
! set the height (zrange) of the surface plot to be related to the x and y extent of the plot  
  max_xy=0d0
  min_xy=1d30
  do ix=1,nx
    do iy=1,ny
      do i=1,2
        max_xy=max(max_xy,abs(xy_values(i,ix,iy)))
        min_xy=min(min_xy,abs(xy_values(i,ix,iy)))
      end do
    end do
  end do
  
  write(*,*)'max_xy=',max_xy
  write(*,*)'min_xy=',min_xy
  
  zrange=(max_xy-min_xy)/4d0
  
  min_value=min(min_Re_value,min_Im_value)
  value_range=max(Re_value_range,Im_value_range)
  
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
  
  do frame=1,n_frames
  
    write(*,*)'Writing frame data',frame
    
! add frame number to base filename

    CALL add_integer_to_filename(filename,frame,filename2)

    frame_filename=trim(filename2)//'.vtk'
    
    write(*,*)'Writing frame to file:',trim(frame_filename)
  
! get the values to plot at this phase point
    do ix=1,nx
      do iy=1,ny
      
        re=dble(cmplx_plot_data(ix,iy))
        im=imag(cmplx_plot_data(ix,iy))
      
        phase=2d0*pi*dble(frame)/dble(n_frames)

        plot_data(ix,iy)=re*cos(phase)-im*sin(phase)
      
      end do ! next iy
    end do ! next ix
 
! Write the data to vtk file

    CALL write_2D_vtk_data(nxplot,nyplot,xy_values,plot_data,min_value,value_range,zrange,frame_filename)    

  end do ! next frame 
  
  DEALLOCATE( xy_values )
  DEALLOCATE( i_values )
  DEALLOCATE( j_values )
  DEALLOCATE( cmplx_plot_data )
  DEALLOCATE( plot_data )
  DEALLOCATE ( data )
  DEALLOCATE ( max_data )
  DEALLOCATE ( min_data )
  DEALLOCATE( column_list )

  RETURN
     
9000 CALL write_line('Error reading data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(filename)
     STOP 
  
END SUBROUTINE animate_2D_complex_data
!
! _______________________________________________________________
!
!    
SUBROUTINE write_2D_vtk_data(nx,ny,xy_values,data,min_value,value_range,zrange,filename)    

  IMPLICIT NONE

  integer	::  nx,ny
  
  real	:: xy_values(1:2,1:nx,1:ny)
  
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
      
      write(20,8000)xy_values(1,ix,iy),xy_values(2,ix,iy),(data(ix,iy)-min_value)*zrange/value_range
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

END SUBROUTINE write_2D_vtk_data
