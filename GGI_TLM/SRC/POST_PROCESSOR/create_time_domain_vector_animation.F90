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
! SUBROUTINE create_time_domain_vector_animation
! SUBROUTINE create_time_domain_near_field_scan
!
! NAME
!    create_time_domain_vector_animation
!
! DESCRIPTION
!    read time domain output over volumes and create a vector animation from a selection of the data 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 16/01/2014 CJS
!
!
SUBROUTINE create_time_domain_vector_animation

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: filename
  
  integer logscale
  parameter (logscale=1)
  
  integer linscale
  parameter (linscale=2)
  
  character*256 output_filename
  character*256 base_filename
  character*256 temp_filename
  character*256 temp_filename2
  character*256 temp_filename3
  character*256 temp_filename4
  character*256 frame_filename
  
  integer n_volumes,volume_loop,volume
  integer output_volume(3)
  
  TYPE(surface_animation_data),allocatable	:: volume_animation(:)
  
  integer			:: n_points,n_quads,n_frames
  integer 			:: output_volume_xyz(3)
  
  real*8 			:: max_data
  real*8 			:: min_data
  real*8 			:: data_range
  
  integer max_point_number

  integer max_quad_number
  
  character*256 :: ipline

  character(len=256)	:: command
  
  logical	:: file_exists
  
  integer 	:: cube_volume,cube_point,face
  integer 	:: ip_point,ip_quad
  integer 	:: point,p,quad
  integer 	:: point1,point2,point3,point4
  real*8 	:: x,y,z
  
  integer 	:: frame
  integer 	:: i_max_data,i_min_data,i_data
  
  integer 	:: i_red,i_green,i_blue
  
  integer i
  character*3 :: scale_string
  integer     :: scale
  real*8 input_real
  real*8 point_data
  
  real*8 :: real_input
  character :: ch
  logical :: scale_data
  
  real*8			:: field(3)
  real*8			:: centre(3)
  real*8			:: cell_size
  real*8			:: vector(3)
  real*8			:: length,diameter
  real*8			:: plot_length
  real*8			:: arrow_points(2,3)
  
  integer 			:: cell,n_cells
  integer			:: n_points_arrow,n_quads_arrow
  integer			:: surface
  integer			:: template_volume
  
  logical	:: found_volume
  integer 	:: number_of_cone_surfaces
  logical	:: scale_vector

! START

  number_of_cone_surfaces=8

! STAGE 1. Open data file

5 write(*,*)
  write(*,*)'volume/Volume field output files:'
  
  command='ls -ltr *.volume_field.tout'
  CALL system(command)

  write(*,*)'Enter the volume field filename'
  read(*,*)filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
   
  base_filename=trim(filename)
    
! STAGE 2. Read volume_field file

!# NUMBER OF OUTPUT volumeS
!           2
!# START volume FIELD FILE TEMPLATE, OUTPUT volume NUMBER:    1
!# Number of points:
!       288
!# Number of faces:
!        72
!# Number of frames:
!        24

  write(*,*)'Reading volume_field file:'
  
  read(local_file_unit,'(A80)'),ipline
  read(local_file_unit,*)n_volumes
  
  ALLOCATE ( volume_animation(1:n_volumes) )
#
! read header information  
  do volume=1,n_volumes
  
    write(*,*)'volume number',volume
    read(local_file_unit,'(A80)'),ipline
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_points
    volume_animation(volume)%n_points=n_points
     
    write(*,*)'Number of points=',volume_animation(volume)%n_points
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_quads
    volume_animation(volume)%n_quads=n_quads
     
    write(*,*)'Number of quads=',volume_animation(volume)%n_quads
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_frames
    volume_animation(volume)%n_frames=n_frames
     
    write(*,*)'Number of frames=',volume_animation(volume)%n_frames

    write(*,*)'Allocating data'
       
    ALLOCATE ( volume_animation(volume)%points(1:n_points,1:3) ) 
    ALLOCATE ( volume_animation(volume)%quads(1:n_quads,1:4) ) 
    ALLOCATE ( volume_animation(volume)%frame_data(1:n_frames,1:n_quads) ) 
    ALLOCATE ( volume_animation(volume)%magnitude_data(1:n_quads) ) 
  
    volume_animation(volume)%points(1:n_points,1:3)=0D0
    volume_animation(volume)%quads(1:n_quads,1:4)=0
    volume_animation(volume)%frame_data(1:n_frames,1:n_quads)=0D0
  
    write(*,*)'Reading points'

    do point=1,n_points
      read(local_file_unit,*)x,y,z
      volume_animation(volume)%points(point,1)=x
      volume_animation(volume)%points(point,2)=y
      volume_animation(volume)%points(point,3)=z
    end do
    
! read quads    
    write(*,*)'Reading quads'

    do quad=1,n_quads
      read(local_file_unit,*)ip_quad,point1,point2,point3,point4
      volume_animation(volume)%quads(quad,1)=point1
      volume_animation(volume)%quads(quad,2)=point2
      volume_animation(volume)%quads(quad,3)=point3
      volume_animation(volume)%quads(quad,4)=point4
    end do
  
    read(local_file_unit,'(A80)'),ipline
  
    write(*,*)'Finished reading fieldsolve geometry data'
    
    volume_animation(volume)%max_data=-1e30
    volume_animation(volume)%min_data=1e30
  
  end do ! read next volume header
 
! STAGE 3. Read fieldsolve frame data

20  CONTINUE
    
!# START OF VOLUME FIELD OUTPUT DATA
!# OUTPUT VOLUME NUMBER
!           1
!# FRAME NUMBER
!           1
!# NUMBER OF FACES
!          72
!# VOLUME FIELD DATA

    read(local_file_unit,'(A80)',end=30),ipline
    
    read(local_file_unit,'(A80)',end=30),ipline
    read(local_file_unit,*,end=30),volume
    
    read(local_file_unit,'(A80)',end=30),ipline
    read(local_file_unit,*,end=30),frame
    
    read(local_file_unit,'(A80)',end=30),ipline
    read(local_file_unit,*,end=30),n_quads
  
    write(*,*)'Reading volume:',volume,' frame:',frame,' n_values:',n_quads
    
    read(local_file_unit,'(A80)',end=30)ipline
    
    do quad=1,n_quads
    
      read(local_file_unit,*,end=30)volume_animation(volume)%frame_data(frame,quad)
      volume_animation(volume)%max_data=	&
              max(volume_animation(volume)%max_data,volume_animation(volume)%frame_data(frame,quad))
      volume_animation(volume)%min_data=	&
              min(volume_animation(volume)%min_data,volume_animation(volume)%frame_data(frame,quad))
	      
    end do ! next quad volume
    
    read(local_file_unit,'(A80)',end=30)ipline

    GOTO 20
  
30 CONTINUE ! jump here when all frame data has been read

100 CONTINUE  ! read surface numbers for the animation
  
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
    write(record_user_inputs_unit,*)output_volume_xyz(i),' Output surface component number',i
    
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
	write(*,*)'Component surfaces have different amounts of data'
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
  write(*,*)'Data range       :',data_range
  
  write(*,*)'Enter the max vector magnitude to plot (or 0 to use the value from the data)'
  read(*,*)real_input
  if (real_input.NE.0d0) then
    data_range=real_input
  end if
  write(record_user_inputs_unit,*)real_input,' max vector magnitude to plot (0 to use the value from the data)'
  
  write(*,*)'Maximum data value:',max_data
  write(*,*)'Minimum data value:',min_data
 
  max_data=data_range

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
  
! there are 24 points associated with each TLM cell
  n_cells=n_points/24   
  
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
    face=0
    do cube_volume=1,n_cells

! get the complex field for this face    
      do i=1,3
      
        surface=output_volume_xyz(i)
        if (surface.ne.0) then
	
	  field(i)=volume_animation(surface)%frame_data(frame,face+1)
	  
        else
	
	  field(i)=0d0
    
        end if
	
      end do
                    
      cell_size =sqrt( ( volume_animation(template_volume)%points(point+1,1)		&
                        -volume_animation(template_volume)%points(point+2,1) )**2+	&
                       ( volume_animation(template_volume)%points(point+1,2)		&
                        -volume_animation(template_volume)%points(point+2,2) )**2+	&
                       ( volume_animation(template_volume)%points(point+1,3)		&
                        -volume_animation(template_volume)%points(point+2,3) )**2 )
      
      do i=1,3
      
        centre(i)=0d0
        do cube_point=1,24
          centre(i)=centre(i)+volume_animation(template_volume)%points(point+cube_point,i)
	end do
	centre(i)=centre(i)/24d0
	
      end do
      
      vector(:)=field(:)
      length=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
      
      if (length.eq.0d0) then
! we must give the plotting routine a non-zero length vector so just send something very small...
 
        vector(1)=max_data/1000d0
	length=max_data/1000d0
 
      end if
      
      if (length.gt.max_data) then
! restrict the length of the 
        vector(:)=vector(:)*max_data/length
        length=max_data
      end if
            
      volume_animation(template_volume)%magnitude_data(cube_volume)=length
      
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
      
      point=point+24
      face=face+6
      
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

  CLOSE(unit=local_file_unit)
    
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
  STOP
9010 continue
  write(*,*)'Error opening frame file'
  write(*,*)'Filename='
  write(*,*)frame_filename
  STOP
  
9020 continue
  write(*,*)'Output volume must be between 0 and ',n_volumes
  STOP

  RETURN
  

  
END SUBROUTINE create_time_domain_vector_animation
!
! NAME
!    create_time_domain_near_field_scan
!
! DESCRIPTION
!    read time domain output over volumes and create a datafile of the form x y z time Field_x Field_y Field_z
!     
! COMMENTS
!    Based on create_time_domain_vector_animation but with a different output format - used initially at least
!    to provide data for time domain near field scanning simulation
!
!    Needs to be tidied up and all the unused bits stripped out...
!
! HISTORY
!
!     started 22/01/2014 CJS
!
!
SUBROUTINE create_time_domain_near_field_scan

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: filename
  
  integer logscale
  parameter (logscale=1)
  
  integer linscale
  parameter (linscale=2)
  
  character*256 output_filename
  character*256 base_filename
  character*256 temp_filename
  character*256 temp_filename2
  character*256 temp_filename3
  character*256 temp_filename4
  character*256 frame_filename
  
  integer n_volumes,volume_loop,volume
  integer output_volume(3)
  
  TYPE(surface_animation_data),allocatable	:: volume_animation(:)
  
  integer			:: n_points,n_quads,n_frames
  integer 			:: output_volume_xyz(3)
  
  real*8 			:: max_data
  real*8 			:: min_data
  real*8 			:: data_range
  
  integer max_point_number

  integer max_quad_number
  
  character*256 :: ipline

  character(len=256)	:: command
  
  logical	:: file_exists
  
  integer 	:: cube_volume,cube_point,face
  integer 	:: ip_point,ip_quad
  integer 	:: point,p,quad
  integer 	:: point1,point2,point3,point4
  real*8 	:: x,y,z
  
  integer 	:: frame
  integer 	:: i_max_data,i_min_data,i_data
  
  integer 	:: i_red,i_green,i_blue
  
  integer i
  character*3 :: scale_string
  integer     :: scale
  real*8 input_real
  real*8 point_data
  
  real*8 :: real_input
  character :: ch
  logical :: scale_data
  
  real*8			:: field(3)
  real*8			:: centre(3)
  real*8			:: cell_size
  real*8			:: vector(3)
  real*8			:: length,diameter
  real*8			:: plot_length
  real*8			:: arrow_points(2,3)
  
  integer 			:: cell,n_cells
  integer			:: n_points_arrow,n_quads_arrow
  integer			:: surface
  integer			:: template_volume
  
  logical	:: found_volume
  integer 	:: number_of_cone_surfaces
  logical	:: scale_vector
  
  real*8	:: output_time,time_zero,time_increment
  
!  logical,parameter	:: test_flag=.TRUE.
  logical,parameter	:: test_flag=.FALSE.

! START

  number_of_cone_surfaces=8

! STAGE 1. Open data file

5 write(*,*)
  write(*,*)'volume/Volume field output files:'
  
  command='ls -ltr *.volume_field.tout'
  CALL system(command)

  write(*,*)'Enter the volume field filename'
  read(*,*)filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
   
  base_filename=trim(filename)
    
! STAGE 2. Read volume_field file

!# NUMBER OF OUTPUT volumeS
!           2
!# START volume FIELD FILE TEMPLATE, OUTPUT volume NUMBER:    1
!# Number of points:
!       288
!# Number of faces:
!        72
!# Number of frames:
!        24

  write(*,*)'Reading volume_field file:'
  
  read(local_file_unit,'(A80)'),ipline
  read(local_file_unit,*)n_volumes
  
  ALLOCATE ( volume_animation(1:n_volumes) )
#
! read header information  
  do volume=1,n_volumes
  
    write(*,*)'volume number',volume
    read(local_file_unit,'(A80)'),ipline
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_points
    volume_animation(volume)%n_points=n_points
     
    write(*,*)'Number of points=',volume_animation(volume)%n_points
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_quads
    volume_animation(volume)%n_quads=n_quads
     
    write(*,*)'Number of quads=',volume_animation(volume)%n_quads
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_frames
    volume_animation(volume)%n_frames=n_frames
     
    write(*,*)'Number of frames=',volume_animation(volume)%n_frames

    write(*,*)'Allocating data'
       
    ALLOCATE ( volume_animation(volume)%points(1:n_points,1:3) ) 
    ALLOCATE ( volume_animation(volume)%quads(1:n_quads,1:4) ) 
    ALLOCATE ( volume_animation(volume)%frame_data(1:n_frames,1:n_quads) ) 
    ALLOCATE ( volume_animation(volume)%magnitude_data(1:n_quads) ) 
  
    volume_animation(volume)%points(1:n_points,1:3)=0D0
    volume_animation(volume)%quads(1:n_quads,1:4)=0
    volume_animation(volume)%frame_data(1:n_frames,1:n_quads)=0D0
  
    write(*,*)'Reading points'

    do point=1,n_points
      read(local_file_unit,*)x,y,z
      volume_animation(volume)%points(point,1)=x
      volume_animation(volume)%points(point,2)=y
      volume_animation(volume)%points(point,3)=z
    end do
    
! read quads    
    write(*,*)'Reading quads'

    do quad=1,n_quads
      read(local_file_unit,*)ip_quad,point1,point2,point3,point4
      volume_animation(volume)%quads(quad,1)=point1
      volume_animation(volume)%quads(quad,2)=point2
      volume_animation(volume)%quads(quad,3)=point3
      volume_animation(volume)%quads(quad,4)=point4
    end do
  
    read(local_file_unit,'(A80)'),ipline
  
    write(*,*)'Finished reading fieldsolve geometry data'
    
    volume_animation(volume)%max_data=-1e30
    volume_animation(volume)%min_data=1e30
  
  end do ! read next volume header
 
! STAGE 3. Read fieldsolve frame data

20  CONTINUE
    
!# START OF VOLUME FIELD OUTPUT DATA
!# OUTPUT VOLUME NUMBER
!           1
!# FRAME NUMBER
!           1
!# NUMBER OF FACES
!          72
!# VOLUME FIELD DATA

    read(local_file_unit,'(A80)',end=30),ipline
    
    read(local_file_unit,'(A80)',end=30),ipline
    read(local_file_unit,*,end=30),volume
    
    read(local_file_unit,'(A80)',end=30),ipline
    read(local_file_unit,*,end=30),frame
    
    read(local_file_unit,'(A80)',end=30),ipline
    read(local_file_unit,*,end=30),n_quads
  
    write(*,*)'Reading volume:',volume,' frame:',frame,' n_values:',n_quads
    
    read(local_file_unit,'(A80)',end=30)ipline
    
    do quad=1,n_quads
    
      read(local_file_unit,*,end=30)volume_animation(volume)%frame_data(frame,quad)
      volume_animation(volume)%max_data=	&
              max(volume_animation(volume)%max_data,volume_animation(volume)%frame_data(frame,quad))
      volume_animation(volume)%min_data=	&
              min(volume_animation(volume)%min_data,volume_animation(volume)%frame_data(frame,quad))
	      
    end do ! next quad volume
    
    read(local_file_unit,'(A80)',end=30)ipline

    GOTO 20
  
30 CONTINUE ! jump here when all frame data has been read

100 CONTINUE  ! read surface numbers for the animation
  
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
    write(record_user_inputs_unit,*)output_volume_xyz(i),' Output surface component number',i
    
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
	write(*,*)'Component surfaces have different amounts of data'
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
  write(*,*)'Data range       :',data_range
  
  write(*,*)'Maximum data value:',max_data
  write(*,*)'Minimum data value:',min_data

  CLOSE(unit=local_file_unit)
  
  write(*,*)'Enter the time of the first output data frame'
  read(*,*)time_zero
  write(record_user_inputs_unit,*)time_zero,' :time of the first output data frame'
  
  write(*,*)'Enter the time interval between output data frames'
  read(*,*)time_increment
  write(record_user_inputs_unit,*)time_increment,' :time interval between output data frames'

  write(*,*)'Enter the filename for the output scan data'
  read(*,'(A256)')output_filename
  write(record_user_inputs_unit,'(A)')trim(output_filename)
  
  n_points=volume_animation(template_volume)%n_points
  n_quads=volume_animation(template_volume)%n_quads
  n_frames=volume_animation(template_volume)%n_frames
  
! there are 24 points associated with each TLM cell
  n_cells=n_points/24   ! this is the number of output points on the scan
    
  write(*,*)'Opening file:',trim(output_filename)
    
  open(UNIT=animation_output_unit,FILE=trim(output_filename),ERR=9010)
  
  write(animation_output_unit,'(A)')	&
'        x                y                z             time           Field_x         Field_y         Field_z' 
    
  do frame=1,n_frames
  
    write(*,*)'Writing frame data',frame
    
    if (test_flag) then
! add frame number to base filename
      CALL add_integer_to_filename(output_filename,frame,frame_filename)
      write(*,*)'Opening file:',trim(frame_filename)   
      open(UNIT=local_file_unit,FILE=trim(frame_filename),ERR=9010)
    end if
    
    output_time=time_zero+(frame-1)*time_increment
	
    point=0
    face=0
    do cube_volume=1,n_cells

! get the vector field for this face    
      do i=1,3
      
        surface=output_volume_xyz(i)
	
        if (surface.ne.0) then	
	  field(i)=volume_animation(surface)%frame_data(frame,face+1)	  
        else	
	  field(i)=0d0    
        end if
	
      end do

! get the coordinate of the centre of the cell      
      do i=1,3
      
        centre(i)=0d0
        do cube_point=1,24
          centre(i)=centre(i)+volume_animation(template_volume)%points(point+cube_point,i)
	end do
	centre(i)=centre(i)/24d0
	
      end do

      write(animation_output_unit,8000)centre(1),centre(2),centre(3),output_time,field(1),field(2),field(3)
      if (test_flag) write(local_file_unit,8000)centre(1),centre(2),centre(3),output_time,field(1),field(2),field(3)
      
8000  format(7E16.6)

! increase the point number by 24 (4 points for each of the 6 cube faces)      
      point=point+24
      
! increase the face number by 6 (6 faces on each cube)    
      face=face+6
    	  
    end do ! next cell
    
    if (test_flag) close(UNIT=local_file_unit)
    
    write(animation_output_unit,*)
    write(animation_output_unit,*)

  end do ! next frame

  close(UNIT=animation_output_unit)

  
! STAGE 7. Close fieldsolve file and deallocate memory
    
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
  STOP
9010 continue
  write(*,*)'Error opening frame file'
  write(*,*)'Filename='
  write(*,*)output_filename
  STOP
  
9020 continue
  write(*,*)'Output volume must be between 0 and ',n_volumes
  STOP

  RETURN
  

  
END SUBROUTINE create_time_domain_near_field_scan
