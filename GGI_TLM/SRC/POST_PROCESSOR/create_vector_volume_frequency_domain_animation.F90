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
! SUBROUTINE create_vector_volume_frequency_domain_animation
!
! NAME
!    create_vector_volume_frequency_domain_animation
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 28/01/2013 CJS
!     bug fix 22/11/2013 CJS n_cells not set correctly when writing animation data to file
!     5/11/2015  CJS allow the process to work on compressed files
!
SUBROUTINE create_vector_volume_frequency_domain_animation

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  character(len=256)	:: filename
  
  integer logscale
  parameter (logscale=1)
  
  integer linscale
  parameter (linscale=2)
  
  character*256 base_filename
  character*256 output_filename
  character*256 temp_filename
  character*256 temp_filename2
  character*256 temp_filename3
  character*256 temp_filename4
  character*256 frame_filename
  
  integer n_volumes,surface_loop,surface
  integer output_volume
  integer template_volume
  integer output_volume_xyz(3)
  
  TYPE(surface_animation_data),allocatable	:: volume_animation(:)
  
  integer			:: n_points,n_quads,n_cells,n_frames
  integer			:: n_points_arrow,n_quads_arrow
  
  real*8 			:: max_data
  real*8 			:: min_data
  real*8 			:: data_range
  
  real*8			:: field(3)
  real*8			:: centre(3)
  real*8			:: cell_size
  real*8			:: vector(3)
  real*8			:: length,diameter
  real*8			:: plot_length
  real*8			:: arrow_points(2,3)
  
  integer max_point_number

  integer max_quad_number
  
  character*256 :: ipline

  character(len=256)	:: command
  
  logical	:: file_exists
  
  integer 	:: ip_point,ip_quad
  integer 	:: point,quad,cell
  integer 	:: point1,point2,point3,point4
  real*8 	:: x,y,z
  
  integer 	:: frame
  integer 	:: i_max_data,i_min_data,i_data
  
  integer 	:: i_red,i_green,i_blue
  
  integer i
  character*3 	:: scale_string
  integer     	:: scale
  real*8     	:: input_real
  real*8     	:: point_data
  
  character 	:: ch
  logical 	:: scale_data
  
  integer 	:: cx,cy,cz
  real*8	:: re,im,mag,phase
  
  logical	:: found_volume
  integer 	:: number_of_cone_surfaces
  logical	:: scale_vector

! compression stuff
  logical	:: compression_flag
  integer	:: len_filename
  character(len=256)	:: filename2
  character*3   :: extn
  
! START

  number_of_cone_surfaces=8

! STAGE 1. Open data file

5 write(*,*)
  write(*,*)'Volume field output files:'
  
  command='ls -ltr *.frequency_output_volume.fout*'
  CALL system(command)

  write(*,*)'Enter the volume field filename'
  read(*,*)filename
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

  if (compression_flag) then
! We must give the filename without the .gz extension here
    filename2=filename(1:len_filename-3)
    CALL open_output_file_read(local_file_unit,filename2,compression_flag)
    base_filename=trim(filename2)
  else
    OPEN(unit=local_file_unit,file=filename)
    base_filename=trim(filename)
  end if
    
! STAGE 2. Read frequency_output_volume file

!# NUMBER OF FREQUENCY OUTPUT VOLUMES:
!2
!# START SURFACE FIELD FILE TEMPLATE, OUTPUT VOLUME NUMBER:    1
!# Number of points:'
!20
!# Number of faces:'
!80
!# Number of frames:'
!24

  write(*,*)'Reading frequency_output_volume file:'
  
  read(local_file_unit,'(A80)'),ipline
  read(local_file_unit,*)n_volumes
  
  ALLOCATE ( volume_animation(1:n_volumes) )
#
! read header information  
  do surface=1,n_volumes
  
    write(*,*)'volume number',surface,' of',n_volumes
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_points
    volume_animation(surface)%n_points=n_points
     
    write(*,*)'Number of points=',volume_animation(surface)%n_points
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_quads
    volume_animation(surface)%n_quads=n_quads
     
    write(*,*)'Number of quads=',volume_animation(surface)%n_quads
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_frames
    volume_animation(surface)%n_frames=n_frames
     
    write(*,*)'Number of frames=',volume_animation(surface)%n_frames
    
    n_cells=n_points/8
    if (n_cells.ne.n_quads/6) then
      write(*,*)'Error in create_vector_volume_frequency_domain_animation'
      write(*,*)'n_cells.ne.n_quads/6'
      write(*,*)'n_points=',n_points,' n_points/8=',n_points/8
      write(*,*)'n_quads =',n_quads, ' n_quads/6 =',n_quads/6
      STOP
    end if

    write(*,*)'Allocating data'
       
    ALLOCATE ( volume_animation(surface)%points(1:n_points,1:3) ) 
    ALLOCATE ( volume_animation(surface)%quads(1:n_quads,1:4) ) 
    ALLOCATE ( volume_animation(surface)%complex_data(1:n_cells) ) 
    ALLOCATE ( volume_animation(surface)%magnitude_data(1:n_cells) ) 
  
    volume_animation(surface)%points(1:n_points,1:3)=0D0
    volume_animation(surface)%quads(1:n_quads,1:4)=0
    volume_animation(surface)%complex_data(1:n_cells)=0D0
    volume_animation(surface)%magnitude_data(1:n_cells)=0D0
  
    write(*,*)'Reading points'

    do point=1,n_points
      read(local_file_unit,*)x,y,z
      volume_animation(surface)%points(point,1)=x
      volume_animation(surface)%points(point,2)=y
      volume_animation(surface)%points(point,3)=z
    end do
    
! read quads    
    write(*,*)'Reading quads'

    do quad=1,n_quads
      read(local_file_unit,*)ip_quad,point1,point2,point3,point4
      volume_animation(surface)%quads(quad,1)=point1
      volume_animation(surface)%quads(quad,2)=point2
      volume_animation(surface)%quads(quad,3)=point3
      volume_animation(surface)%quads(quad,4)=point4
    end do
    
    volume_animation(surface)%max_data=-1e30
    volume_animation(surface)%min_data=1e30
    
! read data    
    write(*,*)'Reading data'
    do quad=1,n_cells
    
      read(local_file_unit,*)cx,cy,cz,re,im,mag
      
      volume_animation(surface)%complex_data(quad)=cmplx(re,im)
      
      volume_animation(surface)%max_data=max(volume_animation(surface)%max_data,mag)
      volume_animation(surface)%min_data=min(volume_animation(surface)%min_data,mag)
      
    end do

    write(*,*)'Maximum data value:',volume_animation(surface)%max_data
    write(*,*)'Minimum data value:',volume_animation(surface)%min_data
    
    write(*,*)'Finished reading fieldsolve geometry data'
      
  end do ! read next volume
 
  write(*,*)' '
  write(*,*)'Number of volumes read:',n_volumes
  write(*,*)' '
  
! read the number of the volume to write to file  

100 CONTINUE  ! read volume numbers for the animation
  
  do i=1,3
  
    if (i.eq.1) then
      write(*,*)'Enter the number of the volume for the field in the x direction or 0 if this component does not exixst'
    else if(i.eq.2) then
      write(*,*)'Enter the number of the volume for the field in the y direction or 0 if this component does not exixst'
    else if(i.eq.3) then
      write(*,*)'Enter the number of the volume for the field in the z direction or 0 if this component does not exixst'
    end if 
    
    read(*,*)output_volume_xyz(i)
    if ( (output_volume_xyz(i).lt.0).OR.(output_volume_xyz(i).gt.n_volumes) )then
      write(*,*)'Output volume must be between 0 and ',n_volumes
      GOTO 100
    end if
    write(record_user_inputs_unit,*)output_volume_xyz(i),' Output volume component number',i
    
  end do

! Check that we have the same amount of data for each field component
! Also work out scaling values for data

  max_data=-1e30
  min_data=1d30

  found_volume=.FALSE.
  do i=1,3
    surface=output_volume_xyz(i)
    if ( (surface.ne.0).AND.(.NOT.found_volume) ) then
! this is the first component of the volume data found    

      template_volume=surface
      n_quads=volume_animation(surface)%n_quads
      n_points=volume_animation(surface)%n_points
      n_frames=volume_animation(surface)%n_frames
      found_volume=.TRUE.
      max_data=max(volume_animation(surface)%max_data,max_data)
      min_data=min(volume_animation(surface)%min_data,min_data)
      
    else if ( (surface.ne.0).AND.(found_volume) ) then
! An output volume has alreeady been found so just check that the data is consistent 
    
      if (  (n_quads.NE.volume_animation(surface)%n_quads).OR.	&
            (n_points.NE.volume_animation(surface)%n_points).OR.	&
            (n_frames.NE.volume_animation(surface)%n_frames ) ) then
	write(*,*)'Error in create_vector_volume_frequency_domain_animation'
	write(*,*)'Component volumes have different amounts of data'
	STOP
      end if
      max_data=max(volume_animation(surface)%max_data,max_data)
      min_data=min(volume_animation(surface)%min_data,min_data)
      
    end if
  end do

  write(*,*)'Maximum data value:',max_data
  write(*,*)'Minimum data value:',min_data
  
  data_range=(max_data-min_data)
  if (data_range.eq.0d0) data_range=1d0

  write(*,*)'Do you want to scale the size of the vector? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  scale_vector=.FALSE.
  if (ch.eq.'y') then
    scale_vector=.TRUE.
  end if

  write(*,*)'Enter the filename for the animation'
  read(*,'(A256)')output_filename
  write(record_user_inputs_unit,'(A)')trim(output_filename)
  
! add surface number to base filename

  temp_filename=trim(output_filename)

  temp_filename3=trim(temp_filename)//'.frame'
  n_points=volume_animation(template_volume)%n_points
  n_quads=volume_animation(template_volume)%n_quads
  n_frames=volume_animation(template_volume)%n_frames
  
  n_cells=n_points/8   ! bug fix 22/11/2013 CJS 
  
  n_quads_arrow=n_cells*number_of_cone_surfaces
  n_points_arrow=n_quads_arrow*4
  
  do frame=1,n_frames
  
    write(*,*)'Writing frame data',frame
    
! add frame number to base filename

    CALL add_integer_to_filename(temp_filename3,frame,temp_filename4)

    frame_filename=trim(temp_filename4)//'.vtk'
    
    write(*,*)'Opening file:',trim(frame_filename)
    
    open(UNIT=animation_output_unit,FILE=trim(frame_filename),ERR=9010)
    
! write header information    
! write vtk header      
    write(animation_output_unit,'(A)')'# vtk DataFile Version 2.0'
    write(animation_output_unit,'(A)')'trim(frame_filename)'
    write(animation_output_unit,'(A)')'ASCII'
    write(animation_output_unit,'(A)')'DATASET POLYDATA'
    write(animation_output_unit,'(A,I10,A)')'POINTS',n_points_arrow,' float'
    
! write cell data 
    point=0
    do cell=1,n_cells

! get the complex field for this face    
      do i=1,3
      
        surface=output_volume_xyz(i)
        if (surface.ne.0) then
    
          re=dble(volume_animation(surface)%complex_data(cell))
          im=dimag(volume_animation(surface)%complex_data(cell))      
          phase=2d0*pi*dble(frame)/dble(n_frames)     
          point_data=re*cos(phase)-im*sin(phase)
	  field(i)=point_data
	  
        else
	
	  field(i)=0d0
    
        end if
	
      end do
      
      do i=1,3
        centre(i)=0.125d0*( volume_animation(template_volume)%points(point+1,i)+	&
                            volume_animation(template_volume)%points(point+2,i)+	&
                            volume_animation(template_volume)%points(point+3,i)+	&
                            volume_animation(template_volume)%points(point+4,i)+	&
                            volume_animation(template_volume)%points(point+5,i)+	&
                            volume_animation(template_volume)%points(point+6,i)+	&
                            volume_animation(template_volume)%points(point+7,i)+	&
                            volume_animation(template_volume)%points(point+8,i) ) 
      end do
                    
      cell_size =sqrt( ( volume_animation(template_volume)%points(point+1,1)		&
                        -volume_animation(template_volume)%points(point+2,1) )**2+	&
                       ( volume_animation(template_volume)%points(point+1,2)		&
                        -volume_animation(template_volume)%points(point+2,2) )**2+	&
                       ( volume_animation(template_volume)%points(point+1,3)		&
                        -volume_animation(template_volume)%points(point+2,3) )**2 )
      
      vector(:)=field(:)
      length=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
      
      if (length.eq.0d0) then
! we must give the plotting routine a non-zero length vector so just send something very small...
 
        vector(1)=max_data/1000d0
	length=max_data/1000d0
 
      end if
            
      volume_animation(template_volume)%magnitude_data(cell)=length
      
      if (scale_vector) then
        vector(:)=vector(:)*cell_size/max_data
        diameter=(cell_size*length/max_data)/8d0
      else
        vector(:)=vector(:)*cell_size/length
        diameter=cell_size/8d0
      end if
       
      arrow_points(1,1:3)=centre(1:3)-vector(1:3)/2d0
      arrow_points(2,1:3)=centre(1:3)+vector(1:3)/2d0
      
      CALL write_cone_points_vtk(arrow_points(1,1),arrow_points(1,2),arrow_points(1,3),	&
                                 arrow_points(2,1),arrow_points(2,2),arrow_points(2,3),	&
			         diameter,number_of_cone_surfaces,animation_output_unit)
      
      point=point+8
      
8000  format(3E14.5)
    	  
    end do ! next cell
  
! write quad data
    write(animation_output_unit,'(A,2I10)')'POLYGONS',n_quads_arrow,n_quads_arrow*5

    point=0
    do quad=1,n_quads_arrow
    
      write(animation_output_unit,8010)4,point,point+1,point+2,point+3
      point=point+4

8010  format(I3,4I8)
      
    end do ! next cell

! write point based data
    write(animation_output_unit,'(A,I10)')'POINT_DATA ',n_points_arrow
    write(animation_output_unit,'(A)')'SCALARS Field_on_cells float 1'
    write(animation_output_unit,'(A)')'LOOKUP_TABLE field_on_cells_table'

    point=0
    do cell=1,n_cells
      
      point_data=volume_animation(template_volume)%magnitude_data(cell)
      
      do i=1,number_of_cone_surfaces*4
        write(animation_output_unit,8020)point_data
      end do
      
8020  format(E14.5)
    	  
    end do ! next quad

    close(UNIT=animation_output_unit)

  end do ! next frame
  
! STAGE 7. Close fieldsolve file and deallocate memory

  if (compression_flag) then
   CALL close_output_file(local_file_unit,filename2,compression_flag)
  else
    CLOSE(unit=local_file_unit)
  end if
    
! Deallocate memory for animation  
    
  do surface=1,n_volumes
    if ( allocated( volume_animation(surface)%points ) ) DEALLOCATE ( volume_animation(surface)%points ) 
    if ( allocated( volume_animation(surface)%quads ) ) DEALLOCATE ( volume_animation(surface)%quads ) 
    if ( allocated( volume_animation(surface)%complex_data ) ) DEALLOCATE ( volume_animation(surface)%complex_data ) 
    if ( allocated( volume_animation(surface)%magnitude_data ) ) DEALLOCATE ( volume_animation(surface)%magnitude_data ) 
  end do
    
  if ( allocated( volume_animation ) ) DEALLOCATE ( volume_animation )  
  
  RETURN
  
9000 continue
  write(*,*)'Error opening output file'
  write(*,*)'Filename='
  write(*,*)base_filename
  stop
9010 continue
  write(*,*)'Error opening frame file'
  write(*,*)'Filename='
  write(*,*)frame_filename
  stop

  RETURN
  

  
END SUBROUTINE create_vector_volume_frequency_domain_animation
